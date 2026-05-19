from __future__ import annotations

"""
OMIX007582 sample-map audit.

This module documents what can be recovered defensibly from the public
Mammal40 methylation release shipped in the repository:

- the beta matrix exposes 643 technical Sentrix barcode/position IDs,
- the public metadata expose 620 biological rows with group/tissue labels,
- no exact join key is present in the public files bundled here.

The audit therefore produces machine-readable inventories and an explicit
summary row rather than inventing a biological sample map by order.
"""

import argparse
import re
import zipfile
from pathlib import Path
from typing import Iterable, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .logging_utils import get_logger

logger = get_logger(__name__)

DEFAULT_METADATA_ID_COLUMNS = (
    "OriginalSampleName",
    "OriginalSampleName.1",
    "sample",
    "sample_id",
)


def _clean_string_set(values: Iterable[object]) -> set[str]:
    cleaned = set()
    for value in values:
        text = str(value).strip()
        if not text or text.lower() in {"nan", "none", "<na>"}:
            continue
        cleaned.add(text)
    return cleaned


def _sample_sheet_like_files(paths: Sequence[str]) -> list[str]:
    keep = []
    for path in paths:
        suffix = Path(path).suffix.lower()
        if suffix in {".csv", ".tsv", ".txt", ".xlsx", ".xls"}:
            keep.append(Path(path).name)
    return sorted(set(keep))


def load_beta_matrix_technical_ids(matrix_path: Path) -> list[str]:
    header = pd.read_csv(matrix_path, nrows=0)
    if len(header.columns) < 2:
        raise ValueError(f"Beta matrix header is malformed: {matrix_path}")
    return [str(col).strip() for col in header.columns[1:]]


def list_idat_prefixes_from_dir(idat_dir: Optional[Path]) -> Tuple[set[str], list[str]]:
    if idat_dir is None or not idat_dir.exists():
        return set(), []

    prefixes = {
        re.sub(r"_(?:Grn|Red)\.idat$", "", path.name, flags=re.IGNORECASE)
        for path in idat_dir.iterdir()
        if path.name.lower().endswith(".idat")
    }
    extras = sorted(
        path.name
        for path in idat_dir.iterdir()
        if not path.name.lower().endswith(".idat")
    )
    return prefixes, extras


def list_idat_prefixes_from_zip(zip_path: Optional[Path]) -> Tuple[set[str], list[str]]:
    if zip_path is None or not zip_path.exists():
        return set(), []

    prefixes: set[str] = set()
    extras: list[str] = []
    with zipfile.ZipFile(zip_path) as handle:
        for name in handle.namelist():
            base = Path(name).name
            if not base:
                continue
            if base.lower().endswith(".idat"):
                prefixes.add(re.sub(r"_(?:Grn|Red)\.idat$", "", base, flags=re.IGNORECASE))
            else:
                extras.append(base)
    return prefixes, sorted(set(extras))


def build_omix007582_sample_map_audit(
    *,
    matrix_path: Path,
    metadata_path: Path,
    idat_dir: Optional[Path] = None,
    zip_path: Optional[Path] = None,
    candidate_cols: Sequence[str] = DEFAULT_METADATA_ID_COLUMNS,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    technical_ids = load_beta_matrix_technical_ids(matrix_path)
    technical_id_set = set(technical_ids)
    metadata = pd.read_csv(metadata_path)
    metadata.columns = [str(col).strip() for col in metadata.columns]

    idat_dir_prefixes, idat_dir_extras = list_idat_prefixes_from_dir(idat_dir)
    zip_prefixes, zip_extras = list_idat_prefixes_from_zip(zip_path)

    overlap_rows = []
    exact_overlap_union: set[str] = set()
    for col in candidate_cols:
        if col not in metadata.columns:
            continue
        values = _clean_string_set(metadata[col].tolist())
        overlap = sorted(technical_id_set & values)
        exact_overlap_union.update(overlap)
        overlap_rows.append(
            {
                "metadata_column": col,
                "n_unique_values": int(len(values)),
                "n_exact_overlaps": int(len(overlap)),
                "example_overlaps": "; ".join(overlap[:5]),
            }
        )

    sample_sheet_like_dir = _sample_sheet_like_files(idat_dir_extras)
    sample_sheet_like_zip = _sample_sheet_like_files(zip_extras)

    if exact_overlap_union:
        mapping_status = "partial_exact_match_available"
        estimable = True
        reason = (
            "At least one exact overlap exists between beta-matrix technical IDs and "
            "metadata identifier columns. Manual review is still required before "
            "promoting any biological sample map."
        )
        evidence_level = 1
    else:
        mapping_status = "technical_ids_only"
        estimable = False
        reason = (
            "No exact overlap exists between the 643 technical beta-matrix IDs and the "
            "620 public metadata rows, and no sample-sheet-like file is present in the "
            "bundled OMIX007582 archive. A defensible biological sample map cannot be "
            "recovered from the current public files without external key metadata."
        )
        evidence_level = 0

    summary = pd.DataFrame(
        [
            {
                "mapping_status": mapping_status,
                "beta_matrix_technical_ids": int(len(technical_ids)),
                "metadata_rows": int(len(metadata)),
                "candidate_columns_checked": "; ".join(
                    [col for col in candidate_cols if col in metadata.columns]
                ),
                "exact_overlap_ids": int(len(exact_overlap_union)),
                "idat_prefixes_in_directory": int(len(idat_dir_prefixes)),
                "idat_prefixes_in_archive": int(len(zip_prefixes)),
                "beta_equals_idat_directory": bool(technical_id_set == idat_dir_prefixes)
                if idat_dir_prefixes
                else False,
                "beta_equals_idat_archive": bool(technical_id_set == zip_prefixes)
                if zip_prefixes
                else False,
                "sample_sheet_like_files_in_directory": int(len(sample_sheet_like_dir)),
                "sample_sheet_like_files_in_archive": int(len(sample_sheet_like_zip)),
                "available": True,
                "estimable": estimable,
                "reason": reason,
                "n_used": int(len(exact_overlap_union)),
                "method": "omix007582_sample_map_audit",
                "ci_low": np.nan,
                "ci_high": np.nan,
                "evidence_level": evidence_level,
            }
        ]
    )

    technical_inventory = pd.DataFrame({"technical_id": technical_ids})
    technical_inventory["sentrix_barcode"] = technical_inventory["technical_id"].str.extract(
        r"^(\d+)_"
    )
    technical_inventory["sentrix_position"] = technical_inventory["technical_id"].str.extract(
        r"_(R\d+C\d+)$"
    )
    technical_inventory["in_beta_matrix"] = True
    technical_inventory["in_idat_directory"] = technical_inventory["technical_id"].isin(idat_dir_prefixes)
    technical_inventory["in_idat_archive"] = technical_inventory["technical_id"].isin(zip_prefixes)
    technical_inventory["mapped_sample_id"] = pd.NA
    technical_inventory["mapping_confidence"] = "none"
    technical_inventory["mapping_status"] = mapping_status
    technical_inventory["reason"] = reason

    metadata_inventory = metadata.copy()
    if "OriginalSampleName" not in metadata_inventory.columns:
        metadata_inventory["OriginalSampleName"] = pd.NA
    if "OriginalSampleName.1" not in metadata_inventory.columns:
        metadata_inventory["OriginalSampleName.1"] = pd.NA
    if "sample" not in metadata_inventory.columns:
        metadata_inventory["sample"] = pd.NA
    metadata_inventory["mapped_technical_id"] = pd.NA
    metadata_inventory["mapping_confidence"] = "none"
    metadata_inventory["mapping_status"] = mapping_status
    metadata_inventory["reason"] = reason

    overlap_audit = pd.DataFrame(overlap_rows)
    if overlap_audit.empty:
        overlap_audit = pd.DataFrame(
            [
                {
                    "metadata_column": "NA",
                    "n_unique_values": 0,
                    "n_exact_overlaps": 0,
                    "example_overlaps": "",
                }
            ]
        )

    return summary, overlap_audit, technical_inventory, metadata_inventory


def write_omix007582_sample_map_audit(
    *,
    matrix_path: Path,
    metadata_path: Path,
    results_dir: Path,
    idat_dir: Optional[Path] = None,
    zip_path: Optional[Path] = None,
) -> dict[str, Path]:
    summary, overlap_audit, technical_inventory, metadata_inventory = build_omix007582_sample_map_audit(
        matrix_path=matrix_path,
        metadata_path=metadata_path,
        idat_dir=idat_dir,
        zip_path=zip_path,
    )

    results_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "summary": results_dir / "omix007582_sample_map_summary.csv",
        "overlap_audit": results_dir / "omix007582_candidate_overlap_audit.csv",
        "technical_inventory": results_dir / "omix007582_technical_id_inventory.csv",
        "metadata_inventory": results_dir / "omix007582_metadata_inventory.csv",
    }
    summary.to_csv(paths["summary"], index=False)
    overlap_audit.to_csv(paths["overlap_audit"], index=False)
    technical_inventory.to_csv(paths["technical_inventory"], index=False)
    metadata_inventory.to_csv(paths["metadata_inventory"], index=False)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit OMIX007582 sample-map recoverability.")
    parser.add_argument(
        "--matrix",
        type=Path,
        default=Path("data/RAW/data/OMIX007582_beta_matrix.csv"),
        help="Path to the OMIX007582 beta matrix with technical IDs.",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=Path("data/RAW/data/OMIX007582-02.csv"),
        help="Path to the OMIX007582 metadata CSV.",
    )
    parser.add_argument(
        "--idat-dir",
        type=Path,
        default=Path("data/RAW/data/OMIX007582_idat"),
        help="Path to the unpacked OMIX007582 IDAT directory.",
    )
    parser.add_argument(
        "--zip-path",
        type=Path,
        default=Path("data/RAW/data/OMIX007582-03.zip"),
        help="Path to the raw OMIX007582 IDAT archive.",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results"),
        help="Directory where audit CSVs will be written.",
    )
    args = parser.parse_args()

    paths = write_omix007582_sample_map_audit(
        matrix_path=args.matrix,
        metadata_path=args.metadata,
        results_dir=args.results_dir,
        idat_dir=args.idat_dir,
        zip_path=args.zip_path,
    )
    for label, path in paths.items():
        logger.info("Wrote %s to %s", label, path)


if __name__ == "__main__":
    main()
