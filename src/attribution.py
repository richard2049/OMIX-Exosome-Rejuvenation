"""
Cross-species attribution and validation helpers for the SRSC pipeline.
"""

from __future__ import annotations

import io
import re
import zipfile
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .clocks import predict_biological_age, train_transcriptomic_clock
from .logging_utils import get_logger

logger = get_logger(__name__)


def _normalize_token(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


def _resolve_tissue_mapping(value: Any, mapping: Optional[Mapping[str, str]] = None) -> str:
    raw = str(value).strip()
    if not raw:
        return ""
    if not mapping:
        return raw
    normalized_map = {_normalize_token(k): v for k, v in mapping.items()}
    return normalized_map.get(_normalize_token(raw), raw)


def parse_mouse_exosome_sample_id(sample_id: str) -> Dict[str, Any]:
    text = str(sample_id).strip()
    match = re.match(
        r"^(?P<age>\d+)_(?P<sex>[A-Za-z])_(?P<tissue>[^_]+)_(?P<arm>[^_]+)_(?P<replicate>\d+)$",
        text,
    )
    if not match:
        raise ValueError(f"Unrecognized OMIX009283 sample ID format: {sample_id}")

    arm_map = {
        "baseline": "Baseline",
        "ctrl": "Ctrl",
        "veh": "Veh",
        "wt": "WT",
        "ges": "GES",
    }
    sex_raw = match.group("sex").upper()
    sex = {"W": "F", "F": "F", "M": "M"}.get(sex_raw, sex_raw)
    arm_raw = match.group("arm")

    return {
        "sample_id": text,
        "age": float(match.group("age")),
        "sex_raw": sex_raw,
        "sex": sex,
        "tissue": match.group("tissue"),
        "arm_raw": arm_raw,
        "arm": arm_map.get(arm_raw.lower(), arm_raw),
        "replicate": int(match.group("replicate")),
    }


def build_mouse_exosome_metadata(sample_ids: Iterable[str]) -> pd.DataFrame:
    rows = []
    for sample_id in sample_ids:
        try:
            parsed = parse_mouse_exosome_sample_id(str(sample_id))
            parsed["parse_ok"] = True
            parsed["reason"] = ""
        except Exception as exc:
            parsed = {
                "sample_id": str(sample_id),
                "age": np.nan,
                "sex_raw": pd.NA,
                "sex": pd.NA,
                "tissue": pd.NA,
                "arm_raw": pd.NA,
                "arm": pd.NA,
                "replicate": np.nan,
                "parse_ok": False,
                "reason": str(exc),
            }
        rows.append(parsed)
    meta = pd.DataFrame(rows)
    if "sample_id" in meta.columns:
        meta["sample_id"] = meta["sample_id"].astype(str)
    return meta


def _effect_with_uncertainty(
    treated: Sequence[float],
    control: Sequence[float],
    *,
    n_bootstrap: int = 1000,
    n_permutations: int = 1000,
    random_state: int = 42,
) -> Dict[str, Any]:
    treated_arr = np.asarray(list(treated), dtype=float)
    control_arr = np.asarray(list(control), dtype=float)
    treated_arr = treated_arr[np.isfinite(treated_arr)]
    control_arr = control_arr[np.isfinite(control_arr)]

    base = {
        "mean_effect": np.nan,
        "ci_low": np.nan,
        "ci_high": np.nan,
        "permutation_p_value": np.nan,
        "n_treated": int(len(treated_arr)),
        "n_control": int(len(control_arr)),
    }
    if len(treated_arr) == 0 or len(control_arr) == 0:
        return base

    rng = np.random.default_rng(random_state)
    observed = float(np.mean(treated_arr) - np.mean(control_arr))
    base["mean_effect"] = observed

    boot = []
    for _ in range(int(n_bootstrap)):
        bt = rng.choice(treated_arr, size=len(treated_arr), replace=True)
        bc = rng.choice(control_arr, size=len(control_arr), replace=True)
        boot.append(float(np.mean(bt) - np.mean(bc)))
    if boot:
        base["ci_low"] = float(np.percentile(boot, 2.5))
        base["ci_high"] = float(np.percentile(boot, 97.5))

    pooled = np.concatenate([treated_arr, control_arr])
    n_t = len(treated_arr)
    perm = []
    for _ in range(int(n_permutations)):
        shuffled = rng.permutation(pooled)
        perm.append(float(np.mean(shuffled[:n_t]) - np.mean(shuffled[n_t:])))
    if perm:
        perm_arr = np.asarray(perm, dtype=float)
        base["permutation_p_value"] = float(np.mean(np.abs(perm_arr) >= abs(observed)))

    return base


def compute_mouse_exosome_tissue_effects(
    expr_log: pd.DataFrame,
    meta: pd.DataFrame,
    *,
    reference_arms: Sequence[str] = ("Baseline",),
    contrasts: Sequence[Tuple[str, str]] = (("GES", "Veh"), ("WT", "Veh"), ("GES", "WT")),
    sample_id_col: str = "sample_id",
    tissue_col: str = "tissue",
    arm_col: str = "arm",
    age_col: str = "age",
    min_reference_samples: int = 20,
    min_reference_ages: int = 4,
    min_samples_per_group: int = 3,
    n_bootstrap: int = 1000,
    n_permutations: int = 1000,
    random_state: int = 42,
) -> pd.DataFrame:
    rows = []

    if expr_log is None or meta is None or expr_log.empty or meta.empty:
        return pd.DataFrame(
            [
                {
                    "tissue": "NA",
                    "contrast": "NA",
                    "treated_arm": "NA",
                    "control_arm": "NA",
                    "available": False,
                    "estimable": False,
                    "reason": "Mouse exosome matrix and/or metadata unavailable.",
                    "n_used": 0,
                    "method": "mouse_tissue_clock_contrast",
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "evidence_level": 0,
                }
            ]
        )

    meta = meta.copy()
    meta[sample_id_col] = meta[sample_id_col].astype(str)

    for tissue, sub in meta.groupby(tissue_col):
        sample_ids = [sid for sid in sub[sample_id_col].tolist() if sid in expr_log.columns]
        if not sample_ids:
            continue

        tissue_expr = expr_log.loc[:, sample_ids]
        sub = sub.set_index(sample_id_col).loc[sample_ids].reset_index()
        sub[age_col] = pd.to_numeric(sub[age_col], errors="coerce")
        ref = sub.loc[sub[arm_col].isin(reference_arms) & sub[age_col].notna()].copy()

        if len(ref) < int(min_reference_samples) or ref[age_col].nunique() < int(min_reference_ages):
            reason = (
                f"Insufficient reference samples for tissue clock in {tissue} "
                f"(n={len(ref)}, unique ages={ref[age_col].nunique()})."
            )
            for treated_arm, control_arm in contrasts:
                rows.append(
                    {
                        "tissue": tissue,
                        "contrast": f"{treated_arm}_vs_{control_arm}",
                        "treated_arm": treated_arm,
                        "control_arm": control_arm,
                        "n_reference": int(len(ref)),
                        "n_reference_ages": int(ref[age_col].nunique()),
                        "n_treated": 0,
                        "n_control": 0,
                        "mean_effect": np.nan,
                        "permutation_p_value": np.nan,
                        "clock_cv_mae": np.nan,
                        "clock_cv_spearman": np.nan,
                        "available": False,
                        "estimable": False,
                        "reason": reason,
                        "n_used": 0,
                        "method": "mouse_tissue_clock_contrast",
                        "ci_low": np.nan,
                        "ci_high": np.nan,
                        "evidence_level": 0,
                    }
                )
            continue

        try:
            clock, _, metrics = train_transcriptomic_clock(
                tissue_expr,
                ref[[sample_id_col, age_col]].copy(),
                age_col=age_col,
                model="ridge",
                n_top_features=None,
                n_splits=min(5, max(2, int(ref[age_col].nunique()))),
                random_state=random_state,
                cv_group_col=None,
            )
            pred = predict_biological_age(clock, tissue_expr, sub[[sample_id_col]].copy())
            pred = pred.rename(columns={"predicted_age": "predicted_age_mouse"})
            sub = sub.merge(pred, on=sample_id_col, how="left")
        except Exception as exc:
            reason = f"Mouse tissue clock failed for {tissue}: {exc}"
            for treated_arm, control_arm in contrasts:
                rows.append(
                    {
                        "tissue": tissue,
                        "contrast": f"{treated_arm}_vs_{control_arm}",
                        "treated_arm": treated_arm,
                        "control_arm": control_arm,
                        "n_reference": int(len(ref)),
                        "n_reference_ages": int(ref[age_col].nunique()),
                        "n_treated": 0,
                        "n_control": 0,
                        "mean_effect": np.nan,
                        "permutation_p_value": np.nan,
                        "clock_cv_mae": np.nan,
                        "clock_cv_spearman": np.nan,
                        "available": False,
                        "estimable": False,
                        "reason": reason,
                        "n_used": 0,
                        "method": "mouse_tissue_clock_contrast",
                        "ci_low": np.nan,
                        "ci_high": np.nan,
                        "evidence_level": 0,
                    }
                )
            continue

        for treated_arm, control_arm in contrasts:
            treated = sub.loc[sub[arm_col] == treated_arm, "predicted_age_mouse"].dropna().tolist()
            control = sub.loc[sub[arm_col] == control_arm, "predicted_age_mouse"].dropna().tolist()
            if len(treated) < int(min_samples_per_group) or len(control) < int(min_samples_per_group):
                rows.append(
                    {
                        "tissue": tissue,
                        "contrast": f"{treated_arm}_vs_{control_arm}",
                        "treated_arm": treated_arm,
                        "control_arm": control_arm,
                        "n_reference": int(len(ref)),
                        "n_reference_ages": int(ref[age_col].nunique()),
                        "n_treated": int(len(treated)),
                        "n_control": int(len(control)),
                        "mean_effect": np.nan,
                        "permutation_p_value": np.nan,
                        "clock_cv_mae": float(metrics.get("MAE", np.nan)),
                        "clock_cv_spearman": float(metrics.get("spearman_r", np.nan)),
                        "available": False,
                        "estimable": False,
                        "reason": (
                            f"Too few samples for {treated_arm} vs {control_arm} in {tissue} "
                            f"(n_treated={len(treated)}, n_control={len(control)})."
                        ),
                        "n_used": int(len(treated) + len(control)),
                        "method": "mouse_tissue_clock_contrast",
                        "ci_low": np.nan,
                        "ci_high": np.nan,
                        "evidence_level": 0,
                    }
                )
                continue

            res = _effect_with_uncertainty(
                treated,
                control,
                n_bootstrap=n_bootstrap,
                n_permutations=n_permutations,
                random_state=random_state,
            )
            rows.append(
                {
                    "tissue": tissue,
                    "contrast": f"{treated_arm}_vs_{control_arm}",
                    "treated_arm": treated_arm,
                    "control_arm": control_arm,
                    "n_reference": int(len(ref)),
                    "n_reference_ages": int(ref[age_col].nunique()),
                    "n_treated": int(res["n_treated"]),
                    "n_control": int(res["n_control"]),
                    "mean_effect": float(res["mean_effect"]),
                    "permutation_p_value": res["permutation_p_value"],
                    "clock_cv_mae": float(metrics.get("MAE", np.nan)),
                    "clock_cv_spearman": float(metrics.get("spearman_r", np.nan)),
                    "available": True,
                    "estimable": True,
                    "reason": "",
                    "n_used": int(res["n_treated"] + res["n_control"]),
                    "method": "mouse_tissue_clock_contrast",
                    "ci_low": res["ci_low"],
                    "ci_high": res["ci_high"],
                    "evidence_level": 1,
                }
            )

    return pd.DataFrame(rows)


def summarize_mouse_exosome_signatures(mouse_effects: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if mouse_effects is None or mouse_effects.empty:
        return pd.DataFrame(
            [
                {
                    "contrast": "NA",
                    "available": False,
                    "estimable": False,
                    "reason": "Mouse exosome effect table unavailable.",
                    "n_used": 0,
                    "method": "mouse_exosome_signature_summary",
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "evidence_level": 0,
                }
            ]
        )

    for contrast, sub in mouse_effects.groupby("contrast"):
        ok = sub.loc[sub["estimable"].astype(bool)].copy()
        if ok.empty:
            rows.append(
                {
                    "contrast": contrast,
                    "n_tissues_estimable": 0,
                    "mean_effect": np.nan,
                    "median_effect": np.nan,
                    "fraction_negative": np.nan,
                    "available": False,
                    "estimable": False,
                    "reason": "No tissues passed mouse exosome estimability filters.",
                    "n_used": 0,
                    "method": "mouse_exosome_signature_summary",
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "evidence_level": 0,
                }
            )
            continue

        effects = ok["mean_effect"].astype(float)
        rows.append(
            {
                "contrast": contrast,
                "n_tissues_estimable": int(len(ok)),
                "mean_effect": float(effects.mean()),
                "median_effect": float(effects.median()),
                "fraction_negative": float((effects < 0).mean()),
                "available": True,
                "estimable": True,
                "reason": "",
                "n_used": int(len(ok)),
                "method": "mouse_exosome_signature_summary",
                "ci_low": float(ok["ci_low"].astype(float).median()),
                "ci_high": float(ok["ci_high"].astype(float).median()),
                "evidence_level": 1,
            }
        )

    return pd.DataFrame(rows)


def _alignment_frame(
    primate_series: pd.Series,
    external_series: pd.Series,
) -> pd.DataFrame:
    df = pd.DataFrame(
        {
            "macaque_effect": pd.to_numeric(primate_series, errors="coerce"),
            "external_effect": pd.to_numeric(external_series, errors="coerce"),
        }
    ).dropna()
    if df.empty:
        return df

    prim_std = float(df["macaque_effect"].std(ddof=0))
    ext_std = float(df["external_effect"].std(ddof=0))
    df["macaque_z"] = (
        (df["macaque_effect"] - df["macaque_effect"].mean()) / prim_std if prim_std > 0 else 0.0
    )
    df["external_z"] = (
        (df["external_effect"] - df["external_effect"].mean()) / ext_std if ext_std > 0 else 0.0
    )
    df["signed_concordance"] = (
        np.sign(df["macaque_effect"]).replace(0, np.nan)
        == np.sign(df["external_effect"]).replace(0, np.nan)
    ).fillna(False).astype(int)
    df["macaque_rank"] = df["macaque_effect"].rank(method="average")
    df["external_rank"] = df["external_effect"].rank(method="average")
    denom = max(1.0, float(len(df) - 1))
    df["rank_concordance"] = 1.0 - (df["macaque_rank"] - df["external_rank"]).abs() / denom
    df["standardized_effect_similarity"] = 1.0 / (
        1.0 + (df["macaque_z"] - df["external_z"]).abs()
    )
    df["residual_component"] = (df["macaque_z"] - df["external_z"]).abs()
    return df


def compute_exosome_alignment_tables(
    primate_effects: pd.DataFrame,
    mouse_effects: pd.DataFrame,
    *,
    tissue_map: Optional[Mapping[str, str]] = None,
    contrasts: Sequence[str] = ("GES_vs_Veh", "WT_vs_Veh"),
    min_common_tissues: int = 3,
    n_bootstrap: int = 1000,
    n_permutations: int = 1000,
    random_state: int = 42,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    by_tissue_rows = []
    summary_rows = []

    if primate_effects is None or primate_effects.empty or mouse_effects is None or mouse_effects.empty:
        stub = pd.DataFrame(
            [
                {
                    "contrast": "NA",
                    "mouse_tissue": "NA",
                    "primate_tissue": "NA",
                    "available": False,
                    "estimable": False,
                    "reason": "Primate and/or mouse effect tables unavailable for exosome alignment.",
                    "n_used": 0,
                    "method": "cross_species_effect_alignment",
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "evidence_level": 0,
                }
            ]
        )
        summary_stub = pd.DataFrame(
            [
                {
                    "contrast": "NA",
                    "n_common_tissues": 0,
                    "available": False,
                    "estimable": False,
                    "reason": "Primate and/or mouse effect tables unavailable for exosome alignment.",
                    "n_used": 0,
                    "method": "cross_species_effect_alignment_summary",
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "evidence_level": 0,
                }
            ]
        )
        return stub, summary_stub

    prim = primate_effects.copy()
    if "tissue" in prim.columns:
        prim = prim.set_index("tissue")
    prim.index = prim.index.map(lambda x: _resolve_tissue_mapping(x, None))
    prim.index.name = "primate_tissue"

    for contrast in contrasts:
        sub = mouse_effects.loc[
            (mouse_effects["contrast"] == contrast)
            & mouse_effects["estimable"].astype(bool)
        ].copy()
        if sub.empty:
            summary_rows.append(
                {
                    "contrast": contrast,
                    "n_common_tissues": 0,
                    "fraction_signed_concordant": np.nan,
                    "mean_rank_concordance": np.nan,
                    "mean_standardized_effect_similarity": np.nan,
                    "residual_component": np.nan,
                    "spearman_rho": np.nan,
                    "pearson_r": np.nan,
                    "preferred_alignment": "",
                    "available": False,
                    "estimable": False,
                    "reason": f"No estimable mouse effects available for contrast {contrast}.",
                    "n_used": 0,
                    "method": "cross_species_effect_alignment_summary",
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "permutation_p_value": np.nan,
                    "evidence_level": 0,
                }
            )
            continue

        sub["primate_tissue"] = sub["tissue"].map(lambda x: _resolve_tissue_mapping(x, tissue_map))
        join = sub.merge(
            prim.reset_index(),
            on="primate_tissue",
            how="inner",
            suffixes=("_mouse", "_macaque"),
        )
        if join.empty:
            join = pd.DataFrame(columns=["tissue", "primate_tissue", "mouse_effect", "macaque_effect"])
        else:
            if "mean_effect_mouse" in join.columns:
                join = join.rename(columns={"mean_effect_mouse": "mouse_effect"})
            elif "mean_effect" in join.columns:
                join["mouse_effect"] = pd.to_numeric(join["mean_effect"], errors="coerce")

            if "mean_effect_macaque" in join.columns:
                join = join.rename(columns={"mean_effect_macaque": "macaque_effect"})
            elif "mean_effect" in prim.columns:
                join["macaque_effect"] = pd.to_numeric(join["mean_effect"], errors="coerce")

        primate_series = pd.Series(
            pd.to_numeric(join.get("macaque_effect", pd.Series(dtype=float)), errors="coerce").values,
            index=join.get("primate_tissue", pd.Series(dtype=str)).astype(str).values,
        )
        mouse_series = pd.Series(
            pd.to_numeric(join.get("mouse_effect", pd.Series(dtype=float)), errors="coerce").values,
            index=join.get("primate_tissue", pd.Series(dtype=str)).astype(str).values,
        )
        aligned = _alignment_frame(primate_series, mouse_series)

        if len(aligned) < int(min_common_tissues):
            reason = (
                f"Only {len(aligned)} common tissues (< {int(min_common_tissues)} required) "
                f"for contrast {contrast}."
            )
            summary_rows.append(
                {
                    "contrast": contrast,
                    "n_common_tissues": int(len(aligned)),
                    "fraction_signed_concordant": np.nan,
                    "mean_rank_concordance": np.nan,
                    "mean_standardized_effect_similarity": np.nan,
                    "residual_component": np.nan,
                    "spearman_rho": np.nan,
                    "pearson_r": np.nan,
                    "preferred_alignment": "",
                    "available": False,
                    "estimable": False,
                    "reason": reason,
                    "n_used": int(len(aligned)),
                    "method": "cross_species_effect_alignment_summary",
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "permutation_p_value": np.nan,
                    "evidence_level": 0,
                }
            )
            by_tissue_rows.append(
                {
                    "contrast": contrast,
                    "mouse_tissue": "NA",
                    "primate_tissue": "NA",
                    "macaque_effect": np.nan,
                    "mouse_effect": np.nan,
                    "signed_concordance": np.nan,
                    "rank_concordance": np.nan,
                    "standardized_effect_similarity": np.nan,
                    "residual_component": np.nan,
                    "global_spearman_rho": np.nan,
                    "global_pearson_r": np.nan,
                    "mouse_permutation_p_value": np.nan,
                    "available": False,
                    "estimable": False,
                    "reason": reason,
                    "n_used": int(len(aligned)),
                    "method": "cross_species_effect_alignment",
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "evidence_level": 0,
                }
            )
            continue

        spearman_rho = float(aligned["macaque_effect"].corr(aligned["external_effect"], method="spearman"))
        pearson_r = float(aligned["macaque_effect"].corr(aligned["external_effect"], method="pearson"))
        observed_similarity = float(aligned["standardized_effect_similarity"].mean())

        rng = np.random.default_rng(random_state)
        boot = []
        arr = aligned[["macaque_effect", "external_effect"]].to_numpy(dtype=float)
        idx = np.arange(len(arr))
        for _ in range(int(n_bootstrap)):
            bidx = rng.choice(idx, size=len(idx), replace=True)
            frame = _alignment_frame(pd.Series(arr[bidx, 0]), pd.Series(arr[bidx, 1]))
            if not frame.empty:
                boot.append(float(frame["standardized_effect_similarity"].mean()))
        ci_low = float(np.percentile(boot, 2.5)) if boot else np.nan
        ci_high = float(np.percentile(boot, 97.5)) if boot else np.nan

        perm = []
        for _ in range(int(n_permutations)):
            shuffled = rng.permutation(arr[:, 1])
            frame = _alignment_frame(pd.Series(arr[:, 0]), pd.Series(shuffled))
            if not frame.empty:
                perm.append(float(frame["standardized_effect_similarity"].mean()))
        perm_p = float(np.mean(np.asarray(perm, dtype=float) >= observed_similarity)) if perm else np.nan

        for tissue_name, row in aligned.iterrows():
            mouse_row = join.loc[join["primate_tissue"] == tissue_name].iloc[0]
            by_tissue_rows.append(
                {
                    "contrast": contrast,
                    "mouse_tissue": str(mouse_row["tissue"]),
                    "primate_tissue": tissue_name,
                    "macaque_effect": float(row["macaque_effect"]),
                    "mouse_effect": float(row["external_effect"]),
                    "signed_concordance": int(row["signed_concordance"]),
                    "rank_concordance": float(row["rank_concordance"]),
                    "standardized_effect_similarity": float(row["standardized_effect_similarity"]),
                    "residual_component": float(row["residual_component"]),
                    "global_spearman_rho": spearman_rho,
                    "global_pearson_r": pearson_r,
                    "mouse_permutation_p_value": float(mouse_row.get("permutation_p_value", np.nan)),
                    "available": True,
                    "estimable": True,
                    "reason": "",
                    "n_used": int(len(aligned)),
                    "method": "cross_species_effect_alignment",
                    "ci_low": float(mouse_row.get("ci_low", np.nan)),
                    "ci_high": float(mouse_row.get("ci_high", np.nan)),
                    "evidence_level": 3,
                }
            )

        summary_rows.append(
            {
                "contrast": contrast,
                "n_common_tissues": int(len(aligned)),
                "fraction_signed_concordant": float(aligned["signed_concordance"].mean()),
                "mean_rank_concordance": float(aligned["rank_concordance"].mean()),
                "mean_standardized_effect_similarity": observed_similarity,
                "residual_component": float(aligned["residual_component"].mean()),
                "spearman_rho": spearman_rho,
                "pearson_r": pearson_r,
                "preferred_alignment": "",
                "available": True,
                "estimable": True,
                "reason": "",
                "n_used": int(len(aligned)),
                "method": "cross_species_effect_alignment_summary",
                "ci_low": ci_low,
                "ci_high": ci_high,
                "permutation_p_value": perm_p,
                "evidence_level": 3,
            }
        )

    summary_df = pd.DataFrame(summary_rows)
    if not summary_df.empty:
        ok = summary_df.loc[summary_df["estimable"].astype(bool)].copy()
        if not ok.empty:
            preferred = str(
                ok.sort_values("mean_standardized_effect_similarity", ascending=False).iloc[0]["contrast"]
            )
            summary_df["preferred_alignment"] = preferred

    return pd.DataFrame(by_tissue_rows), summary_df


def assign_group_stage_age(
    groups: pd.Series,
    group_age_map: Mapping[str, float],
) -> pd.Series:
    normalized_map = {_normalize_token(k): float(v) for k, v in group_age_map.items()}

    def _map_one(value: Any) -> float:
        key = _normalize_token(value)
        if key in normalized_map:
            return normalized_map[key]
        if key.startswith("y"):
            return normalized_map.get("yc", np.nan)
        if key.startswith("m"):
            return normalized_map.get("mc", np.nan)
        if key.startswith("o"):
            return normalized_map.get("ov", np.nan)
        return np.nan

    return groups.map(_map_one).astype(float)


def summarize_multimodal_concordance(
    transcript_effects: pd.DataFrame,
    methylation_effects: pd.DataFrame,
    *,
    min_common_tissues: int = 3,
    n_bootstrap: int = 1000,
    n_permutations: int = 1000,
    random_state: int = 42,
) -> pd.DataFrame:
    if transcript_effects is None or transcript_effects.empty or methylation_effects is None or methylation_effects.empty:
        return pd.DataFrame(
            [
                {
                    "available": False,
                    "estimable": False,
                    "reason": "Transcriptomic and/or methylation rejuvenation tables unavailable.",
                    "n_used": 0,
                    "method": "multimodal_concordance",
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "evidence_level": 0,
                }
            ]
        )

    tx = transcript_effects.copy()
    me = methylation_effects.copy()
    if "tissue" in tx.columns:
        tx = tx.set_index("tissue")
    if "tissue" in me.columns:
        me = me.set_index("tissue")

    frame = _alignment_frame(
        tx["effect_median"] if "effect_median" in tx.columns else tx["mean_effect"],
        me["effect_median"] if "effect_median" in me.columns else me["mean_effect"],
    )
    if len(frame) < int(min_common_tissues):
        return pd.DataFrame(
            [
                {
                    "n_common_tissues": int(len(frame)),
                    "fraction_signed_concordant": np.nan,
                    "mean_rank_concordance": np.nan,
                    "mean_standardized_effect_similarity": np.nan,
                    "residual_component": np.nan,
                    "spearman_rho": np.nan,
                    "pearson_r": np.nan,
                    "available": False,
                    "estimable": False,
                    "reason": (
                        f"Only {len(frame)} common tissues (< {int(min_common_tissues)} required) "
                        "for multimodal concordance."
                    ),
                    "n_used": int(len(frame)),
                    "method": "multimodal_concordance",
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "permutation_p_value": np.nan,
                    "evidence_level": 0,
                }
            ]
        )

    observed_similarity = float(frame["standardized_effect_similarity"].mean())
    rng = np.random.default_rng(random_state)
    boot = []
    arr = frame[["macaque_effect", "external_effect"]].to_numpy(dtype=float)
    idx = np.arange(len(arr))
    for _ in range(int(n_bootstrap)):
        bidx = rng.choice(idx, size=len(idx), replace=True)
        boot_frame = _alignment_frame(pd.Series(arr[bidx, 0]), pd.Series(arr[bidx, 1]))
        if not boot_frame.empty:
            boot.append(float(boot_frame["standardized_effect_similarity"].mean()))
    ci_low = float(np.percentile(boot, 2.5)) if boot else np.nan
    ci_high = float(np.percentile(boot, 97.5)) if boot else np.nan

    perm = []
    for _ in range(int(n_permutations)):
        shuffled = rng.permutation(arr[:, 1])
        perm_frame = _alignment_frame(pd.Series(arr[:, 0]), pd.Series(shuffled))
        if not perm_frame.empty:
            perm.append(float(perm_frame["standardized_effect_similarity"].mean()))
    perm_p = float(np.mean(np.asarray(perm, dtype=float) >= observed_similarity)) if perm else np.nan

    return pd.DataFrame(
        [
            {
                "n_common_tissues": int(len(frame)),
                "common_tissues": ";".join(frame.index.astype(str)),
                "fraction_signed_concordant": float(frame["signed_concordance"].mean()),
                "mean_rank_concordance": float(frame["rank_concordance"].mean()),
                "mean_standardized_effect_similarity": observed_similarity,
                "residual_component": float(frame["residual_component"].mean()),
                "spearman_rho": float(frame["macaque_effect"].corr(frame["external_effect"], method="spearman")),
                "pearson_r": float(frame["macaque_effect"].corr(frame["external_effect"], method="pearson")),
                "available": True,
                "estimable": True,
                "reason": "",
                "n_used": int(len(frame)),
                "method": "multimodal_concordance",
                "ci_low": ci_low,
                "ci_high": ci_high,
                "permutation_p_value": perm_p,
                "evidence_level": 1,
            }
        ]
    )


def read_subset_sample_info(zip_path: Path) -> pd.DataFrame:
    with zipfile.ZipFile(zip_path) as zf:
        members = [name for name in zf.namelist() if name.endswith("sample.info.csv")]
        if not members:
            raise FileNotFoundError(f"No sample.info.csv found in {zip_path.name}")
        with zf.open(members[0]) as handle:
            text = io.TextIOWrapper(handle, encoding="utf-8", newline="")
            return pd.read_csv(text)


def build_subset_validation_table(
    sample_info: pd.DataFrame,
    *,
    subset_name: str,
    bulk_tissue: str,
    group_map: Mapping[str, str],
    bulk_effects: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    df = sample_info.copy()
    raw_group = df["group"].astype(str).str.strip() if "group" in df.columns else pd.Series([], dtype=str)
    repo_group = raw_group.map(lambda x: group_map.get(x, x))
    counts = repo_group.value_counts().to_dict()

    bulk_row = None
    if bulk_effects is not None and not bulk_effects.empty:
        tmp = bulk_effects.copy()
        if "tissue" in tmp.columns:
            match = tmp.loc[tmp["tissue"].astype(str) == str(bulk_tissue)]
            if not match.empty:
                bulk_row = match.iloc[0]
        elif bulk_tissue in tmp.index:
            bulk_row = tmp.loc[bulk_tissue]

    row = {
        "subset": subset_name,
        "bulk_tissue": bulk_tissue,
        "n_samples_total": int(len(df)),
        "n_y_c": int(counts.get("Y_C", 0)),
        "n_m_c": int(counts.get("M_C", 0)),
        "n_o_c": int(counts.get("O_C", 0)),
        "n_o_v": int(counts.get("O_V", 0)),
        "n_o_wt": int(counts.get("O_WT", 0)),
        "n_o_ges": int(counts.get("O_GES", 0)),
        "group_distinction_preserved": bool((counts.get("O_C", 0) > 0) and (counts.get("O_V", 0) > 0)),
        "intervention_groups_present": bool(
            counts.get("O_V", 0) > 0 and counts.get("O_WT", 0) > 0 and counts.get("O_GES", 0) > 0
        ),
        "bulk_effect_median": np.nan,
        "bulk_ci_low": np.nan,
        "bulk_ci_high": np.nan,
        "validation_status": "metadata_only",
        "available": True,
        "estimable": False,
        "reason": (
            "Subset sample-sheet validation integrated. Expression-level pseudobulk validation "
            "for zipped single-cell matrices is not yet implemented."
        ),
        "n_used": int(len(df)),
        "method": "subset_sample_info_audit",
        "ci_low": np.nan,
        "ci_high": np.nan,
        "evidence_level": 1,
    }
    if bulk_row is not None:
        row["bulk_effect_median"] = float(
            bulk_row["effect_median"] if "effect_median" in bulk_row else bulk_row.get("mean_effect", np.nan)
        )
        row["bulk_ci_low"] = float(bulk_row.get("ci_low", np.nan))
        row["bulk_ci_high"] = float(bulk_row.get("ci_high", np.nan))
    else:
        row["reason"] = (
            "Subset sample-sheet validation integrated, but no matching bulk tissue effect "
            "was available for direct comparison."
        )
        row["evidence_level"] = 0

    return pd.DataFrame([row])
