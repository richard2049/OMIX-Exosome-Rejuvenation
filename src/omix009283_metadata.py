from __future__ import annotations

"""
Build a deterministic OMIX009283 sample-metadata table from matrix headers.

The released OMIX009283 matrix encodes age, sex, tissue, intervention arm,
and replicate index in the sample IDs. This helper reconstructs that metadata
into a reusable CSV instead of relying on one-off local scripts.
"""

import argparse
from pathlib import Path

import pandas as pd

from .attribution import build_mouse_exosome_metadata
from .logging_utils import get_logger

logger = get_logger(__name__)

FEATURE_ID_COLUMNS = {"gene", "genes", "gene_id", "symbol"}


def extract_omix009283_sample_ids(matrix_path: Path) -> list[str]:
    header = pd.read_csv(matrix_path, sep="\t", nrows=0)
    sample_ids = [
        str(col).strip()
        for col in header.columns
        if str(col).strip().lower() not in FEATURE_ID_COLUMNS
    ]
    if not sample_ids:
        raise ValueError(f"No sample IDs detected in {matrix_path}")
    return sample_ids


def build_omix009283_metadata_table(matrix_path: Path) -> pd.DataFrame:
    sample_ids = extract_omix009283_sample_ids(matrix_path)
    metadata = build_mouse_exosome_metadata(sample_ids)
    metadata = metadata.sort_values("sample_id").reset_index(drop=True)
    return metadata


def write_omix009283_metadata(matrix_path: Path, output_path: Path) -> Path:
    metadata = build_omix009283_metadata_table(matrix_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(output_path, index=False)

    parse_ok = int(metadata["parse_ok"].fillna(False).sum()) if "parse_ok" in metadata.columns else 0
    logger.info(
        "Wrote OMIX009283 metadata to %s (%d rows, %d parsed successfully).",
        output_path,
        len(metadata),
        parse_ok,
    )
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Build deterministic OMIX009283 metadata from matrix headers.")
    parser.add_argument(
        "--matrix",
        type=Path,
        default=Path("data/RAW/data/OMIX009283-01.txt"),
        help="Path to the released OMIX009283 matrix.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/RAW/data/OMIX009283_metadata.csv"),
        help="Output CSV path for the rebuilt metadata table.",
    )
    args = parser.parse_args()
    write_omix009283_metadata(args.matrix, args.output)


if __name__ == "__main__":
    main()
