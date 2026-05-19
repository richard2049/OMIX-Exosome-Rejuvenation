"""
rejuvenation.py

Helpers to compute rejuvenation-related metrics from the clock output.

This module provides:
  - delta_age = predicted_age - chronological_age
  - global rejuvenation effect (treated vs control, with bootstrap CI)
  - tissue-level rejuvenation summary (per-tissue effect sizes)
  - simple expression-based tissue effects (mean treated-control shift)

It assumes:
  - Sample-level metadata includes group labels, tissues, and ages
  - Clock predictions are already computed and stored in the metadata
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple, Sequence

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
from .logging_utils import get_logger

logger = get_logger(__name__)


def _benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    p = pd.to_numeric(p_values, errors="coerce").to_numpy(dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    mask = np.isfinite(p)
    if mask.sum() == 0:
        return pd.Series(out, index=p_values.index, dtype=float)

    pv = p[mask]
    order = np.argsort(pv)
    ranked = pv[order]
    m = float(len(ranked))
    q = ranked * m / (np.arange(1, len(ranked) + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)

    out_idx = np.where(mask)[0][order]
    out[out_idx] = q
    return pd.Series(out, index=p_values.index, dtype=float)


def compute_delta_age(
    meta: pd.DataFrame,
    pred_age_col: str,
    chrono_age_col: str,
    out_col: str = "delta_age",
) -> pd.DataFrame:
    """
    Add a delta_age column = predicted_age - chronological_age.

    Fails cleanly if required columns are missing or too many values are NaN.
    """
    if pred_age_col not in meta.columns:
        raise ValueError(f"Missing predicted age column: {pred_age_col}")
    if chrono_age_col not in meta.columns:
        raise ValueError(f"Missing chronological age column: {chrono_age_col}")

    df = meta.copy()
    df[out_col] = df[pred_age_col] - df[chrono_age_col]

    # Basic NaN sanity check
    valid_frac = df[out_col].notna().mean()
    if valid_frac < 0.5:
        raise ValueError(
            f"Too few valid delta_age values ({valid_frac:.2%}). "
            "Check age prediction and metadata."
        )
    return df


def _group_effect(
    df: pd.DataFrame,
    group_col: str,
    value_col: str,
    control_labels: List[str],
    treated_labels: List[str],
    min_per_group: int = 4,
    n_bootstrap: int = 2000,
    random_state: int = 42,
) -> Optional[Dict]:
    """
    Overall or tissue-level rejuvenation effect:
    difference in medians (treated - control) and bootstrap CI.
    Returns None if there is no minimum power.
    """
    if group_col not in df.columns:
        return None

    g = df.dropna(subset=[group_col, value_col])

    is_control = g[group_col].isin(control_labels)
    is_treated = g[group_col].isin(treated_labels)

    ctrl = g.loc[is_control, value_col].values
    trt = g.loc[is_treated, value_col].values

    if len(ctrl) < min_per_group or len(trt) < min_per_group:
        return None

    rng = np.random.default_rng(random_state)
    diffs = []
    for _ in range(n_bootstrap):
        s_ctrl = rng.choice(ctrl, size=len(ctrl), replace=True)
        s_trt = rng.choice(trt, size=len(trt), replace=True)
        diffs.append(np.median(s_trt) - np.median(s_ctrl))

    diffs = np.array(diffs)
    effect = float(np.median(diffs))
    ci_low, ci_high = np.percentile(diffs, [2.5, 97.5])

    return {
        "n_ctrl": int(len(ctrl)),
        "n_trt": int(len(trt)),
        "effect_median": effect,
        "ci_low": float(ci_low),
        "ci_high": float(ci_high),
    }


def summarize_rejuvenation_by_tissue(
    meta_with_delta: pd.DataFrame,
    tissue_col: str,
    group_col: str,
    value_col: str = "delta_age",
    control_labels: List[str] = ("Control", "WTC", "Saline"),
    treated_labels: List[str] = ("SRC", "V", "GES"),
    min_per_group: int = 4,
    n_bootstrap: int = 2000,
    random_state: int = 42,
) -> pd.DataFrame:
    """
    Returns a rejuvenation by tissue table.
    Filters automatically tissues without effect.
    """
    if tissue_col not in meta_with_delta.columns:
        raise ValueError(f"Missing tissue column: {tissue_col}")
    if group_col not in meta_with_delta.columns:
        raise ValueError(f"Missing group column: {group_col}")

    rows = []
    for tissue, sub in meta_with_delta.groupby(tissue_col):
        eff = _group_effect(
            sub,
            group_col=group_col,
            value_col=value_col,
            control_labels=list(control_labels),
            treated_labels=list(treated_labels),
            min_per_group=min_per_group,
            n_bootstrap=n_bootstrap,
            random_state=random_state,
        )
        if eff is None:
            continue
        eff_row = {
            "tissue": tissue,
            **eff,
        }
        rows.append(eff_row)

    return pd.DataFrame(rows)


def summarise_global_rejuvenation(
    meta_with_delta: pd.DataFrame,
    group_col: str,
    value_col: str = "delta_age",
    control_labels: List[str] = ("Control", "WTC", "Saline"),
    treated_labels: List[str] = ("SRC", "V", "GES"),
    min_per_group: int = 6,
    n_bootstrap: int = 2000,
    random_state: int = 42,
) -> Optional[Dict]:
    """
    Global effect size (not stratified by tissue).

    Returns None if the minimum power requirements are not met.
    """
    return _group_effect(
        meta_with_delta,
        group_col=group_col,
        value_col=value_col,
        control_labels=list(control_labels),
        treated_labels=list(treated_labels),
        min_per_group=min_per_group,
        n_bootstrap=n_bootstrap,
        random_state=random_state,
    )

def summarize_tissue_expression_effects(
    expr: pd.DataFrame,
    meta: pd.DataFrame,
    tissue_col: str,
    group_col: str,
    control_labels: Sequence[str],
    treated_labels: Sequence[str],
    outcome_col: str = "delta_age",
    min_per_group: int = 2,
    covariate_cols: Optional[Sequence[str]] = ("age", "sex", "batch"),
) -> pd.DataFrame:
    """
    Summarize tissue-wise treatment effects on a sample-level outcome using
    per-tissue linear models:

        outcome ~ treated + covariates

    This avoids pseudo-replication from treating genes as independent units
    for tissue-level inference.

    Parameters
    ----------
    expr : DataFrame
        Gene expression matrix (genes x samples). Used for sample alignment and
        gene-count metadata; inference is performed at sample level.
    meta : DataFrame
        Sample metadata. Must contain tissue_col, group_col, sample_id and outcome_col.
    tissue_col : str
        Column name for tissue/organ.
    group_col : str
        Column name for experimental group (e.g. Y_C, O_C, O_GES...).
    control_labels : list-like
        Labels considered as control.
    treated_labels : list-like
        Labels considered as treated.
    min_per_group : int
        Minimum number of samples per group within a tissue to compute effects.

    Returns one row per tissue with covariate-adjusted treatment-effect summaries.
    """

    out_cols = [
        "tissue",
        "n_ctrl",
        "n_trt",
        "n_samples",
        "n_used",
        "n_genes_modeled",
        "method",
        "covariates_used",
        "top_genes",
        "mean_effect",
        "median_effect",
        "effect_se",
        "ci_low",
        "ci_high",
        "p_value",
        "fdr_q_value",
        "available",
        "estimable",
        "reason",
    ]

    required_cols = {tissue_col, group_col, "sample_id", outcome_col}
    missing = required_cols - set(meta.columns)
    if missing:
        logger.warning(
            "summarize_tissue_expression_effects: missing required "
            "metadata columns: %s. Returning empty DataFrame.",
            ", ".join(sorted(missing)),
        )
        return pd.DataFrame(columns=out_cols)

    # Drop rows without tissue/group/outcome
    meta = meta.copy()
    meta = meta.loc[
        meta[tissue_col].notna()
        & meta[group_col].notna()
        & pd.to_numeric(meta[outcome_col], errors="coerce").notna()
    ]
    if meta.empty:
        logger.warning(
            "summarize_tissue_expression_effects: no rows with tissue/group/outcome. Returning empty."
        )
        return pd.DataFrame(columns=out_cols)

    # Keep only control + treated labels
    valid_labels = list(control_labels) + list(treated_labels)
    meta = meta.loc[meta[group_col].isin(valid_labels)]
    if meta.empty:
        logger.warning(
            "summarize_tissue_expression_effects: no rows with group in %s. Returning empty.",
            valid_labels,
        )
        return pd.DataFrame(columns=out_cols)

    # Use sample_id as index to align with expr columns
    meta["sample_id"] = meta["sample_id"].astype(str)
    meta = meta.set_index("sample_id")

    # Align expression columns to metadata
    common_ids = expr.columns.intersection(meta.index)
    if len(common_ids) < (2 * min_per_group):
        logger.warning(
            "summarize_tissue_expression_effects: only %d common samples between expr and meta. Returning empty.",
            len(common_ids),
        )
        return pd.DataFrame(columns=out_cols)

    expr = expr.loc[:, common_ids]
    meta = meta.loc[common_ids]

    def _encode_covariate(name: str, s: pd.Series) -> Optional[pd.DataFrame]:
        s_num = pd.to_numeric(s, errors="coerce")
        if s_num.notna().mean() >= 0.8 and s_num.nunique(dropna=True) > 1:
            return pd.DataFrame({name: s_num.astype(float)})

        name_low = name.lower()
        if ("sex" in name_low) or ("gender" in name_low):
            s_low = s.astype("string").str.strip().str.lower()
            mapped = s_low.map(
                {
                    "m": 1.0,
                    "male": 1.0,
                    "f": 0.0,
                    "female": 0.0,
                }
            )
            if mapped.notna().mean() >= 0.8 and mapped.nunique(dropna=True) > 1:
                return pd.DataFrame({name: mapped.astype(float)})

        # Generic categorical encoding (e.g., batch)
        s_cat = s.astype("string").str.strip().replace({"": pd.NA, "nan": pd.NA, "None": pd.NA, "<NA>": pd.NA})
        if s_cat.notna().mean() < 0.8 or s_cat.nunique(dropna=True) <= 1:
            return None
        if s_cat.nunique(dropna=True) > 12:
            return None
        dummies = pd.get_dummies(s_cat, prefix=name, drop_first=True, dtype=float)
        if dummies.shape[1] == 0:
            return None
        return dummies

    rows = []
    for tissue, sub_meta in meta.groupby(tissue_col):
        # Boolean masks in this tissue
        is_ctrl = sub_meta[group_col].isin(control_labels)
        is_trt = sub_meta[group_col].isin(treated_labels)

        n_ctrl = int(is_ctrl.sum())
        n_trt = int(is_trt.sum())

        if n_ctrl < min_per_group or n_trt < min_per_group:
            continue

        ids = [s for s in sub_meta.index.tolist() if s in expr.columns]
        if len(ids) < (2 * min_per_group):
            continue
        sub = sub_meta.loc[ids].copy()
        sub["treated_binary"] = sub[group_col].isin(treated_labels).astype(float)
        sub[outcome_col] = pd.to_numeric(sub[outcome_col], errors="coerce")

        design = pd.DataFrame(index=sub.index)
        design["intercept"] = 1.0
        design["treated_binary"] = sub["treated_binary"].astype(float)
        used_covariates: List[str] = []

        for cov in list(covariate_cols or []):
            if cov not in sub.columns:
                continue
            encoded = _encode_covariate(cov, sub[cov])
            if encoded is None:
                continue
            for cname in encoded.columns:
                design[cname] = encoded[cname].astype(float)
            used_covariates.append(cov)

        design["__y__"] = sub[outcome_col].astype(float)
        design = design.dropna(axis=0, how="any")
        if design.empty:
            continue
        if design["treated_binary"].nunique() < 2:
            continue

        keep_ids = design.index
        if len(keep_ids) < (2 * min_per_group):
            continue

        n_ctrl_model = int((design["treated_binary"] == 0).sum())
        n_trt_model = int((design["treated_binary"] == 1).sum())
        if n_ctrl_model < min_per_group or n_trt_model < min_per_group:
            continue

        X_df = design.drop(columns=["__y__"])
        y = design["__y__"].to_numpy(dtype=float)
        X = X_df.to_numpy(dtype=float)
        n, p = X.shape
        if n <= (p + 1):
            continue

        beta = np.linalg.lstsq(X, y, rcond=None)[0]
        y_hat = X @ beta
        resid = y - y_hat
        dof = n - p
        if dof <= 0:
            continue

        sigma2 = float(np.dot(resid, resid) / dof)
        xtx_inv = np.linalg.pinv(X.T @ X)
        var_treat = float(max(sigma2 * xtx_inv[1, 1], 0.0))
        effect = float(beta[1])
        effect_se = float(np.sqrt(var_treat)) if np.isfinite(var_treat) else np.nan

        if np.isfinite(effect_se) and effect_se > 0:
            t_stat = effect / effect_se
            p_value = float(2.0 * stats.t.sf(np.abs(t_stat), dof))
            t_crit = float(stats.t.ppf(0.975, dof))
            ci_low = float(effect - t_crit * effect_se)
            ci_high = float(effect + t_crit * effect_se)
        else:
            p_value = np.nan
            ci_low, ci_high = np.nan, np.nan

        rows.append(
            {
                "tissue": tissue,
                "n_ctrl": n_ctrl_model,
                "n_trt": n_trt_model,
                "n_samples": int(n),
                "n_used": int(n),
                "n_genes_modeled": int(expr.shape[0]),
                "method": "ols_sample_level_outcome_by_tissue",
                "covariates_used": ",".join(used_covariates) if used_covariates else "none",
                "top_genes": "",
                "mean_effect": effect,
                "median_effect": effect,
                "effect_se": effect_se,
                "ci_low": float(ci_low) if not np.isnan(ci_low) else np.nan,
                "ci_high": float(ci_high) if not np.isnan(ci_high) else np.nan,
                "p_value": p_value,
                "fdr_q_value": np.nan,
                "available": True,
                "estimable": True,
                "reason": "",
            }
        )

    if not rows:
        logger.warning(
            "summarize_tissue_expression_effects: no tissues passed the filters (min_per_group=%d). "
            "Returning empty DataFrame.",
            min_per_group,
        )
        return pd.DataFrame(columns=out_cols)

    out = pd.DataFrame(rows, columns=out_cols)
    out["fdr_q_value"] = _benjamini_hochberg(out["p_value"])
    return out

def compute_plasma_biomarkers(
    plasma_expr: pd.DataFrame,
    plasma_meta: pd.DataFrame,
    outcome: pd.Series,
    min_non_nan_frac: float = 0.7,
    min_pairs: int = 8,
    n_bootstrap: int = 200,
    random_state: int = 42,
    min_sign_agreement: float = 0.8,
    stability_top_k: int = 500,
) -> pd.DataFrame:
    """
    Rank plasma proteins by association with an outcome (e.g. rejuvenation score).

    Parameters
    ----------
    plasma_expr : DataFrame
        Rows = proteins, cols = samples (matching plasma_meta index or sample_id).
    plasma_meta : DataFrame
        Must be index-aligned or contain a column to align with plasma_expr columns.
    outcome : Series
        Numeric phenotype per sample (e.g. rejuvenation score, state index, etc.).
        Index MUST be sample IDs matching plasma_expr columns.
    min_non_nan_frac : float
        Minimum fraction of non-NaN values per protein to include in the analysis.

    Returns
    -------
    DataFrame with association and stability columns.
    """
    # Ensuring alignment
    common = plasma_expr.columns.intersection(outcome.index)
    if len(common) < 3:
        raise ValueError(f"Too few samples with both plasma and outcome: {len(common)}")

    expr = plasma_expr.loc[:, common]
    y = outcome.loc[common].astype(float)

    # Filtered by missingness
    non_nan_frac = expr.notna().mean(axis=1)
    expr = expr.loc[non_nan_frac >= min_non_nan_frac]
    if expr.empty:
        raise ValueError("All plasma features dropped due to missingness.")

    records = []
    y_values = y.values.astype(float)
    for row_id, (protein, row) in enumerate(expr.iterrows()):
        x = row.values.astype(float)
        mask = ~np.isnan(x) & ~np.isnan(y_values)
        n_pairs_eff = int(mask.sum())
        if n_pairs_eff < int(min_pairs):
            continue
        xv = x[mask]
        yv = y_values[mask]
        r, p = stats.spearmanr(xv, yv)
        if np.isnan(r):
            continue
        records.append((row_id, protein, r, p, n_pairs_eff))

    if not records:
        return pd.DataFrame(
            columns=[
                "protein",
                "spearman_r",
                "pval",
                "qval",
                "abs_r",
                "direction",
                "n_pairs",
                "rho_ci_low",
                "rho_ci_high",
                "sign_agreement",
                "stable_association",
                "stability_tested",
            ]
        )

    df = pd.DataFrame(
        records,
        columns=[
            "__row_id",
            "protein",
            "spearman_r",
            "pval",
            "n_pairs",
        ],
    )

    # Simple FDR (Benjamini-Hochberg)
    df = df.sort_values("pval").reset_index(drop=True)
    m = len(df)
    df["qval"] = df["pval"] * m / (df.index + 1)
    df["qval"] = df["qval"].clip(upper=1.0)

    df["abs_r"] = df["spearman_r"].abs()
    df["direction"] = np.where(df["spearman_r"] > 0, "pro-aging", "pro-rejuvenation")
    df["rho_ci_low"] = np.nan
    df["rho_ci_high"] = np.nan
    df["sign_agreement"] = np.nan
    df["stable_association"] = False
    df["stability_tested"] = False

    # Bootstrap stability is computed on top-ranked proteins to keep runtime bounded.
    n_test = int(min(max(0, int(stability_top_k)), len(df)))
    if int(n_bootstrap) > 0 and n_test > 0:
        rng = np.random.default_rng(random_state)
        top_idx = df.sort_values("abs_r", ascending=False).head(n_test).index.tolist()
        y_full = y.astype(float)
        for idx in top_idx:
            row_id = int(df.at[idx, "__row_id"])
            if row_id < 0 or row_id >= expr.shape[0]:
                continue
            x = expr.iloc[row_id, :].astype(float).values
            mask = ~np.isnan(x) & ~np.isnan(y_full.values)
            n_pairs_eff = int(mask.sum())
            if n_pairs_eff < int(min_pairs):
                continue

            xv = x[mask]
            yv = y_full.values[mask]
            bs_r = []
            for _ in range(int(n_bootstrap)):
                bi = rng.integers(0, n_pairs_eff, size=n_pairs_eff)
                xb = xv[bi]
                yb = yv[bi]
                if np.nanstd(xb) == 0 or np.nanstd(yb) == 0:
                    continue
                rb = stats.spearmanr(xb, yb).correlation
                if np.isfinite(rb):
                    bs_r.append(float(rb))
            if not bs_r:
                continue

            ci_low = float(np.percentile(bs_r, 2.5))
            ci_high = float(np.percentile(bs_r, 97.5))
            pos_frac = float(np.mean(np.asarray(bs_r) > 0))
            neg_frac = float(np.mean(np.asarray(bs_r) < 0))
            sign_agreement = float(max(pos_frac, neg_frac))
            ci_excludes_zero = bool(ci_low > 0 or ci_high < 0)
            stable_association = bool(
                ci_excludes_zero
                and np.isfinite(sign_agreement)
                and sign_agreement >= float(min_sign_agreement)
            )

            df.at[idx, "rho_ci_low"] = ci_low
            df.at[idx, "rho_ci_high"] = ci_high
            df.at[idx, "sign_agreement"] = sign_agreement
            df.at[idx, "stable_association"] = stable_association
            df.at[idx, "stability_tested"] = True
    return df.drop(columns=["__row_id"], errors="ignore")

def simple_mediation_bootstrap(
    df: pd.DataFrame,
    x_col: str,
    m_col: str,
    y_col: str,
    n_boot: int = 1000,
    seed: int = 42,
) -> dict:
    """
    Simple (non-parametric) mediation bootstrap:
    X -> M -> Y, with Y also adjusted by X.

    Returns:
      - total_effect
      - direct_effect
      - indirect_effect
      - bootstrap 95% CIs
    """
    df = df[[x_col, m_col, y_col]].dropna().copy()
    if df.shape[0] < 10:
        logger.warning(
            "simple_mediation_bootstrap: too few samples (n=%d). Returning NaNs.",
            df.shape[0],
        )
        return {
            "n_samples": int(df.shape[0]),
            "total_effect": np.nan,
            "direct_effect": np.nan,
            "indirect_effect": np.nan,
            "total_ci_low": np.nan,
            "total_ci_high": np.nan,
            "direct_ci_low": np.nan,
            "direct_ci_high": np.nan,
            "indirect_ci_low": np.nan,
            "indirect_ci_high": np.nan,
        }

    rng = np.random.default_rng(seed)

    X = df[[x_col]].to_numpy(dtype=float)
    M = df[[m_col]].to_numpy(dtype=float)
    Y = df[[y_col]].to_numpy(dtype=float)

    # a) Total effect: Y ~ X
    reg_total = LinearRegression().fit(X, Y)
    total_effect = float(reg_total.coef_[0, 0])

    # b) Direct and mediated effects:
    #    1) M ~ X  -> coefficient a
    reg_a = LinearRegression().fit(X, M)
    a = float(reg_a.coef_[0, 0])

    #    2) Y ~ X + M -> coefficient b (M), coefficient c' (X)
    XM = np.concatenate([X, M], axis=1)
    reg_b = LinearRegression().fit(XM, Y)
    b = float(reg_b.coef_[0, 1])  # coefficient for M
    direct_effect = float(reg_b.coef_[0, 0])  # coefficient for X
    indirect_effect = a * b

    total_samples = []
    direct_samples = []
    indirect_samples = []

    for _ in range(n_boot):
        idx = rng.integers(0, df.shape[0], size=df.shape[0])
        Xb = X[idx]
        Mb = M[idx]
        Yb = Y[idx]

        reg_total_b = LinearRegression().fit(Xb, Yb)
        total_b = float(reg_total_b.coef_[0, 0])

        reg_a_b = LinearRegression().fit(Xb, Mb)
        a_b = float(reg_a_b.coef_[0, 0])

        XM_b = np.concatenate([Xb, Mb], axis=1)
        reg_b_b = LinearRegression().fit(XM_b, Yb)
        b_b = float(reg_b_b.coef_[0, 1])
        direct_b = float(reg_b_b.coef_[0, 0])
        indirect_b = a_b * b_b

        total_samples.append(total_b)
        direct_samples.append(direct_b)
        indirect_samples.append(indirect_b)

    def ci(arr):
        lo, hi = np.percentile(arr, [2.5, 97.5])
        return float(lo), float(hi)

    t_lo, t_hi = ci(total_samples)
    d_lo, d_hi = ci(direct_samples)
    i_lo, i_hi = ci(indirect_samples)

    out = {
        "n_samples": int(df.shape[0]),
        "total_effect": total_effect,
        "direct_effect": direct_effect,
        "indirect_effect": indirect_effect,
        "total_ci_low": t_lo,
        "total_ci_high": t_hi,
        "direct_ci_low": d_lo,
        "direct_ci_high": d_hi,
        "indirect_ci_low": i_lo,
        "indirect_ci_high": i_hi,
    }
    logger.info("Mediation bootstrap done: %s", out)
    return out

def run_simple_causal_model(
    prim_meta: pd.DataFrame,
    prim_plasma_meta: pd.DataFrame,
    cfg,
) -> Optional[pd.DataFrame]:
    """
    Phase 4 helper: build a simple causal model:
    group_binary -> plasma_state_score -> rejuvenation_score
    on matched animals.
    """
    if "animal_id" not in prim_meta.columns or "animal_id" not in prim_plasma_meta.columns:
        logger.warning("run_simple_causal_model: animal_id missing; cannot align plasma and tissue.")
        return None

    if "plasma_state_score" not in prim_plasma_meta.columns:
        logger.warning("run_simple_causal_model: plasma_state_score missing; skipping mediation.")
        return None

    if "rejuvenation_score" not in prim_meta.columns:
        logger.warning("run_simple_causal_model: rejuvenation_score missing; skipping mediation.")
        return None

    meta = prim_meta.copy()
    if "group_binary" not in meta.columns:
        meta["group_binary"] = (meta["group"] == cfg.primate_treated_label).astype(int)

    merged = meta.merge(
        prim_plasma_meta[["animal_id", "plasma_state_score"]],
        on="animal_id",
        how="inner",
    )
    merged = merged.dropna(subset=["group_binary", "plasma_state_score", "rejuvenation_score"])

    if merged.shape[0] < cfg.min_samples_for_mediation:
        logger.warning(
            "run_simple_causal_model: too few matched samples for mediation (n=%d).",
            merged.shape[0],
        )
        return None

    med = simple_mediation_bootstrap(
        merged,
        x_col="group_binary",
        m_col="plasma_state_score",
        y_col="rejuvenation_score",
        n_boot=cfg.mediation_bootstrap,
        seed=cfg.random_seed,
    )

    return pd.DataFrame([med])
