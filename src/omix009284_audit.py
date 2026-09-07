from __future__ import annotations

"""
Read-only OMIX009284 audit.

The current plan explicitly defers OMIX009284 integration until the repository
has a non-destructive structural audit. This module characterizes the released
files without attempting to fold them into the main attribution pipeline.

Key questions answered here:
- What assay/export type do the raw files resemble?
- Can file-level naming be linked to sample/group labels defensibly?
- Does the dataset materially strengthen the current cross-tissue attribution
  question, or is it better treated as a future PBMC-focused follow-up?
"""

import argparse
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .logging_utils import get_logger
from .reason_codes import (
    OMIX009284_PBMC_ONLY_FOR_ATTRIBUTION,
    AUTHOR_KEY_EXOSOME_UPTAKE,
    annotate_reason_fields,
)

logger = get_logger(__name__)

BARCODE_SUFFIX_RE = re.compile(r"^[A-Z]+-\d+_(\d+)$", re.IGNORECASE)
MATRIX_ROLE = "gene_by_cell_matrix"
METADATA_ROLE = "cell_metadata_table"
UNKNOWN_ROLE = "unknown"
REQUIRED_METADATA_COLUMNS = {"orig.ident", "group", "sample", "tissue", "cell"}


def _read_preview_lines(path: Path) -> tuple[list[str], list[str]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        first_line = handle.readline().rstrip("\r\n")
        second_line = handle.readline().rstrip("\r\n")
    first_tokens = first_line.split("\t") if first_line else []
    second_tokens = second_line.split("\t") if second_line else []
    return first_tokens, second_tokens


def _looks_like_matrix_header(tokens: list[str]) -> bool:
    if len(tokens) < 2:
        return False
    preview = tokens[: min(25, len(tokens))]
    if not preview:
        return False
    return all(BARCODE_SUFFIX_RE.match(token) for token in preview)


def _extract_matrix_suffix(tokens: Iterable[str]) -> str:
    suffixes = {
        match.group(1)
        for token in tokens
        if (match := BARCODE_SUFFIX_RE.match(str(token).strip()))
    }
    if len(suffixes) == 1:
        return next(iter(suffixes))
    return ""


def _classify_file(path: Path) -> dict[str, object]:
    header_tokens, second_tokens = _read_preview_lines(path)
    role = UNKNOWN_ROLE
    assay_hint = "undetermined"
    barcode_suffix = ""

    if REQUIRED_METADATA_COLUMNS.issubset(set(header_tokens)):
        role = METADATA_ROLE
        assay_hint = "seurat_style_single_cell_metadata"
    elif _looks_like_matrix_header(header_tokens):
        role = MATRIX_ROLE
        assay_hint = "seurat_style_single_cell_expression_matrix"
        barcode_suffix = _extract_matrix_suffix(header_tokens)

    return {
        "file_name": path.name,
        "file_size_bytes": int(path.stat().st_size),
        "header_fields": int(len(header_tokens)),
        "first_token": header_tokens[0] if header_tokens else "",
        "second_row_first_token": second_tokens[0] if second_tokens else "",
        "role": role,
        "assay_hint": assay_hint,
        "barcode_suffix": barcode_suffix,
        "n_cells_header": int(len(header_tokens)) if role == MATRIX_ROLE else np.nan,
    }


def _load_metadata_table(metadata_path: Path) -> pd.DataFrame:
    metadata = pd.read_csv(metadata_path, sep="\t", index_col=0)
    metadata.columns = [str(col).strip() for col in metadata.columns]
    missing = REQUIRED_METADATA_COLUMNS - set(metadata.columns)
    if missing:
        raise ValueError(
            f"Metadata table is missing required OMIX009284 columns: {sorted(missing)}"
        )
    return metadata


def _unique_join(values: pd.Series) -> str:
    return ";".join(sorted({str(value) for value in values.dropna().astype(str)}))


def _build_suffix_sample_map(
    metadata: pd.DataFrame,
    matrix_inventory: pd.DataFrame,
) -> pd.DataFrame:
    barcode_values = metadata["cell"].astype(str)
    suffix = barcode_values.str.extract(r"_(\d+)$")[0]
    mapping = (
        pd.DataFrame(
            {
                "suffix": suffix,
                "sample": metadata["sample"].astype(str),
                "group": metadata["group"].astype(str),
                "tissue": metadata["tissue"].astype(str),
            }
        )
        .dropna(subset=["suffix"])
        .groupby("suffix", as_index=False)
        .agg(
            n_cells_metadata=("sample", "size"),
            n_samples=("sample", "nunique"),
            sample=("sample", _unique_join),
            n_groups=("group", "nunique"),
            group=("group", _unique_join),
            n_tissues=("tissue", "nunique"),
            tissue=("tissue", _unique_join),
        )
    )

    matrix_lookup = matrix_inventory[["file_name", "barcode_suffix", "n_cells_header"]].rename(
        columns={"barcode_suffix": "suffix"}
    )
    mapping = mapping.merge(matrix_lookup, on="suffix", how="left")
    mapping["has_expression_matrix"] = mapping["file_name"].notna()
    mapping["cell_count_match"] = mapping["n_cells_metadata"] == mapping["n_cells_header"]
    mapping["mapping_confidence"] = np.where(
        (mapping["n_samples"] == 1)
        & (mapping["n_groups"] == 1)
        & (mapping["n_tissues"] == 1)
        & mapping["has_expression_matrix"]
        & mapping["cell_count_match"],
        "high",
        "ambiguous",
    )
    mapping["mapping_status"] = np.where(
        mapping["mapping_confidence"] == "high",
        "recoverable_from_metadata_suffix",
        "requires_manual_review",
    )
    mapping = mapping.sort_values("suffix", key=lambda s: s.astype(int)).reset_index(drop=True)
    return mapping


def build_omix009284_audit(
    *,
    data_dir: Path,
    pattern: str = "OMIX009284-*.txt",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    paths = sorted(data_dir.glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No OMIX009284 files matched {pattern} under {data_dir}")

    file_inventory = pd.DataFrame([_classify_file(path) for path in paths])
    matrix_inventory = file_inventory[file_inventory["role"] == MATRIX_ROLE].copy()
    metadata_rows = file_inventory[file_inventory["role"] == METADATA_ROLE].copy()

    metadata_summary = {
        "metadata_rows": 0,
        "unique_samples": 0,
        "unique_groups": 0,
        "unique_tissues": 0,
        "group_labels": "",
        "tissue_labels": "",
    }
    sample_mapping = pd.DataFrame(
        columns=[
            "suffix",
            "n_cells_metadata",
            "n_samples",
            "sample",
            "n_groups",
            "group",
            "n_tissues",
            "tissue",
            "file_name",
            "n_cells_header",
            "has_expression_matrix",
            "cell_count_match",
            "mapping_confidence",
            "mapping_status",
        ]
    )

    if len(metadata_rows) == 1:
        metadata_path = data_dir / str(metadata_rows.iloc[0]["file_name"])
        metadata = _load_metadata_table(metadata_path)
        metadata_summary = {
            "metadata_rows": int(len(metadata)),
            "unique_samples": int(metadata["sample"].nunique(dropna=True)),
            "unique_groups": int(metadata["group"].nunique(dropna=True)),
            "unique_tissues": int(metadata["tissue"].nunique(dropna=True)),
            "group_labels": _unique_join(metadata["group"]),
            "tissue_labels": _unique_join(metadata["tissue"]),
        }
        sample_mapping = _build_suffix_sample_map(metadata, matrix_inventory)

    matrix_suffixes = {
        str(value)
        for value in matrix_inventory["barcode_suffix"].dropna().astype(str)
        if str(value).strip()
    }
    mapping_suffixes = {
        str(value)
        for value in sample_mapping["suffix"].dropna().astype(str)
        if str(value).strip()
    }
    exact_suffix_coverage = bool(matrix_suffixes) and matrix_suffixes == mapping_suffixes
    exact_mapping = bool(
        exact_suffix_coverage
        and not sample_mapping.empty
        and sample_mapping["mapping_confidence"].eq("high").all()
    )

    if exact_mapping:
        mapping_status = "recoverable_from_metadata_suffix"
    elif metadata_rows.empty:
        mapping_status = "metadata_table_missing"
    else:
        mapping_status = "partial_or_ambiguous"

    reason = (
        "Read-only audit indicates that OMIX009284 is a PBMC single-cell RNA export "
        "with 32 per-sample gene-by-cell matrices plus one Seurat-style cell metadata "
        "table. Sample suffixes can be linked back to sample/group labels through the "
        "metadata table, but the dataset remains PBMC-only and therefore does not "
        "materially strengthen the current cross-tissue exosome-alignment question."
    )

    summary = annotate_reason_fields(pd.DataFrame(
        [
            {
                "dataset_scope": "PBMC single-cell RNA export",
                "n_files": int(len(file_inventory)),
                "n_expression_matrices": int((file_inventory["role"] == MATRIX_ROLE).sum()),
                "n_metadata_tables": int((file_inventory["role"] == METADATA_ROLE).sum()),
                "n_unknown_files": int((file_inventory["role"] == UNKNOWN_ROLE).sum()),
                "metadata_rows": int(metadata_summary["metadata_rows"]),
                "unique_samples": int(metadata_summary["unique_samples"]),
                "unique_groups": int(metadata_summary["unique_groups"]),
                "unique_tissues": int(metadata_summary["unique_tissues"]),
                "group_labels": str(metadata_summary["group_labels"]),
                "tissue_labels": str(metadata_summary["tissue_labels"]),
                "suffix_mapping_status": mapping_status,
                "matrix_suffixes_match_metadata": exact_suffix_coverage,
                "material_for_current_attribution": False,
                "available": True,
                "estimable": False,
                "reason": reason,
                "reason_code": OMIX009284_PBMC_ONLY_FOR_ATTRIBUTION,
                "missing_author_key": AUTHOR_KEY_EXOSOME_UPTAKE,
                "n_used": int(len(file_inventory)),
                "method": "omix009284_read_only_audit",
                "ci_low": np.nan,
                "ci_high": np.nan,
                "evidence_level": 0,
            }
        ]
    ))

    return summary, file_inventory, sample_mapping


def write_omix009284_audit(
    *,
    data_dir: Path,
    results_dir: Path,
    pattern: str = "OMIX009284-*.txt",
) -> dict[str, Path]:
    summary, file_inventory, sample_mapping = build_omix009284_audit(
        data_dir=data_dir,
        pattern=pattern,
    )
    results_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "summary": results_dir / "omix009284_audit_summary.csv",
        "file_inventory": results_dir / "omix009284_file_inventory.csv",
        "sample_mapping": results_dir / "omix009284_suffix_sample_map.csv",
    }
    summary.to_csv(paths["summary"], index=False)
    file_inventory.to_csv(paths["file_inventory"], index=False)
    sample_mapping.to_csv(paths["sample_mapping"], index=False)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description="Run a read-only structural audit of OMIX009284.")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("data/RAW/data"),
        help="Directory containing the OMIX009284 raw text exports.",
    )
    parser.add_argument(
        "--pattern",
        type=str,
        default="OMIX009284-*.txt",
        help="Filename glob used to locate OMIX009284 files under data-dir.",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results"),
        help="Directory where audit CSVs will be written.",
    )
    args = parser.parse_args()

    paths = write_omix009284_audit(
        data_dir=args.data_dir,
        results_dir=args.results_dir,
        pattern=args.pattern,
    )
    for label, path in paths.items():
        logger.info("Wrote %s to %s", label, path)


if __name__ == "__main__":
    main()
