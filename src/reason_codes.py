from __future__ import annotations

"""
Stable reason codes for interpretation-facing pipeline outputs.

The human-readable ``reason`` column remains the primary explanation for a row.
These codes add machine-readable triage so non-estimable results can be traced
back to either a concrete missing author-side key or a non-metadata limitation.
"""

from typing import Any, Mapping, Optional

import pandas as pd

OK = "OK"
CONFIG_DISABLED = "CONFIG_DISABLED"
DATA_FILE_UNAVAILABLE = "DATA_FILE_UNAVAILABLE"
SCHEMA_REQUIREMENT_MISSING = "SCHEMA_REQUIREMENT_MISSING"
LOW_SAMPLE_SIZE_CAUTION = "LOW_SAMPLE_SIZE_CAUTION"
INSUFFICIENT_GROUP_SAMPLE_SIZE = "INSUFFICIENT_GROUP_SAMPLE_SIZE"
INSUFFICIENT_SHARED_TISSUE_OVERLAP = "INSUFFICIENT_SHARED_TISSUE_OVERLAP"
INSUFFICIENT_LINKED_ANIMALS = "INSUFFICIENT_LINKED_ANIMALS"
PLASMA_BULK_ANIMAL_LINKAGE_MISSING = "PLASMA_BULK_ANIMAL_LINKAGE_MISSING"
PLASMA_LINKAGE_CONFIDENCE_MISSING = "PLASMA_LINKAGE_CONFIDENCE_MISSING"
PLASMA_ANIMAL_LINKAGE_COLLISION = "PLASMA_ANIMAL_LINKAGE_COLLISION"
EXOSOME_CARGO_RECIPIENT_LINKAGE_MISSING = "EXOSOME_CARGO_RECIPIENT_LINKAGE_MISSING"
OMIX007582_SENTRIX_SAMPLE_SHEET_MISSING = "OMIX007582_SENTRIX_SAMPLE_SHEET_MISSING"
OMIX009284_PBMC_ONLY_FOR_ATTRIBUTION = "OMIX009284_PBMC_ONLY_FOR_ATTRIBUTION"
CONTROL_STAGE_DESIGN_KEY_REQUIRED = "CONTROL_STAGE_DESIGN_KEY_REQUIRED"
QC_EXCLUSION_TABLE_REQUIRED = "QC_EXCLUSION_TABLE_REQUIRED"
SRSC_EXOSOME_PREP_METADATA_REQUIRED = "SRSC_EXOSOME_PREP_METADATA_REQUIRED"
SUBSET_PSEUDOBULK_NOT_IMPLEMENTED = "SUBSET_PSEUDOBULK_NOT_IMPLEMENTED"
NO_ESTIMABLE_FEATURES = "NO_ESTIMABLE_FEATURES"
METHOD_RUNTIME_FAILURE = "METHOD_RUNTIME_FAILURE"
NON_ESTIMABLE_UNSPECIFIED = "NON_ESTIMABLE_UNSPECIFIED"

AUTHOR_KEY_NOT_REQUIRED = "not_required"
AUTHOR_KEY_REVIEW_REQUIRED = "review_required"
AUTHOR_KEY_PLASMA_TO_ANIMAL = (
    "omix007581_plasma_sample_id_to_animal_id_group_sex_age"
)
AUTHOR_KEY_BULK_SAMPLE_MAP = (
    "omix007580_bulk_rnaseq_sample_id_to_animal_id_tissue_group_sex_age"
)
AUTHOR_KEY_PLASMA_BULK_LINKAGE = (
    f"{AUTHOR_KEY_PLASMA_TO_ANIMAL};{AUTHOR_KEY_BULK_SAMPLE_MAP}"
)
AUTHOR_KEY_EXOSOME_CARGO_RECIPIENT = (
    "exosome_preparation_donor_cargo_treatment_recipient_outcome_linkage"
)
AUTHOR_KEY_OMIX007582_SENTRIX = (
    "omix007582_sentrix_barcode_position_to_biological_sample_tissue_group_sex_age"
)
AUTHOR_KEY_QC_EXCLUSION = "rna_plasma_methylation_single_cell_qc_exclusion_tables"
AUTHOR_KEY_SRSC_PREP = (
    "srsc_exosome_source_foxo3a_validation_dose_schedule_batch_recipient_metadata"
)
AUTHOR_KEY_EXOSOME_UPTAKE = (
    "direct_exosome_uptake_biodistribution_or_cargo_transfer_recipient_tissue_key"
)
AUTHOR_KEY_CONTROL_STAGE = "treatment_era_vs_stage_control_design_key"

MISSING_AUTHOR_KEY_BY_REASON_CODE = {
    OK: "",
    CONFIG_DISABLED: AUTHOR_KEY_NOT_REQUIRED,
    DATA_FILE_UNAVAILABLE: AUTHOR_KEY_NOT_REQUIRED,
    SCHEMA_REQUIREMENT_MISSING: AUTHOR_KEY_NOT_REQUIRED,
    LOW_SAMPLE_SIZE_CAUTION: AUTHOR_KEY_NOT_REQUIRED,
    INSUFFICIENT_GROUP_SAMPLE_SIZE: AUTHOR_KEY_NOT_REQUIRED,
    INSUFFICIENT_SHARED_TISSUE_OVERLAP: AUTHOR_KEY_NOT_REQUIRED,
    INSUFFICIENT_LINKED_ANIMALS: AUTHOR_KEY_PLASMA_BULK_LINKAGE,
    PLASMA_BULK_ANIMAL_LINKAGE_MISSING: AUTHOR_KEY_PLASMA_BULK_LINKAGE,
    PLASMA_LINKAGE_CONFIDENCE_MISSING: AUTHOR_KEY_PLASMA_TO_ANIMAL,
    PLASMA_ANIMAL_LINKAGE_COLLISION: AUTHOR_KEY_PLASMA_BULK_LINKAGE,
    EXOSOME_CARGO_RECIPIENT_LINKAGE_MISSING: AUTHOR_KEY_EXOSOME_CARGO_RECIPIENT,
    OMIX007582_SENTRIX_SAMPLE_SHEET_MISSING: AUTHOR_KEY_OMIX007582_SENTRIX,
    OMIX009284_PBMC_ONLY_FOR_ATTRIBUTION: AUTHOR_KEY_EXOSOME_UPTAKE,
    CONTROL_STAGE_DESIGN_KEY_REQUIRED: AUTHOR_KEY_CONTROL_STAGE,
    QC_EXCLUSION_TABLE_REQUIRED: AUTHOR_KEY_QC_EXCLUSION,
    SRSC_EXOSOME_PREP_METADATA_REQUIRED: AUTHOR_KEY_SRSC_PREP,
    SUBSET_PSEUDOBULK_NOT_IMPLEMENTED: AUTHOR_KEY_NOT_REQUIRED,
    NO_ESTIMABLE_FEATURES: AUTHOR_KEY_NOT_REQUIRED,
    METHOD_RUNTIME_FAILURE: AUTHOR_KEY_REVIEW_REQUIRED,
    NON_ESTIMABLE_UNSPECIFIED: AUTHOR_KEY_REVIEW_REQUIRED,
}


def _clean_optional(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass
    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "<na>"}:
        return ""
    return text


def _as_bool(value: Any, default: bool = True) -> bool:
    if value is None:
        return default
    try:
        if pd.isna(value):
            return default
    except TypeError:
        pass
    if isinstance(value, str):
        text = value.strip().lower()
        if text in {"true", "1", "yes", "y"}:
            return True
        if text in {"false", "0", "no", "n"}:
            return False
    return bool(value)


def infer_reason_code(
    *,
    reason: str = "",
    method: str = "",
    estimable: bool = True,
) -> str:
    reason_l = str(reason or "").lower()
    method_l = str(method or "").lower()
    text = f"{reason_l} {method_l}"

    if estimable and not reason_l:
        return OK
    if "disabled by config" in text or "block disabled" in text:
        return CONFIG_DISABLED
    if "omix007582" in text or "mammal40" in text or "sentrix" in text:
        return OMIX007582_SENTRIX_SAMPLE_SHEET_MISSING
    if "methylation_stage" in text or "age_proxy_years" in text:
        return OMIX007582_SENTRIX_SAMPLE_SHEET_MISSING
    if "technical beta-matrix" in text or "biological sample map" in text:
        return OMIX007582_SENTRIX_SAMPLE_SHEET_MISSING
    if "pbmc-only" in text or "pbmc single-cell" in text or "pbmc focused" in text:
        return OMIX009284_PBMC_ONLY_FOR_ATTRIBUTION
    if "exosome cargo" in text or "donor" in text and "recipient" in text:
        return EXOSOME_CARGO_RECIPIENT_LINKAGE_MISSING
    if "uptake" in text or "biodistribution" in text or "cargo-transfer" in text:
        return OMIX009284_PBMC_ONLY_FOR_ATTRIBUTION
    if "foxo3a" in text or "dose" in text or "injection schedule" in text:
        return SRSC_EXOSOME_PREP_METADATA_REQUIRED
    if "qc/exclusion" in text or "exclusion table" in text:
        return QC_EXCLUSION_TABLE_REQUIRED
    if "stage controls" in text or "treatment-era controls" in text:
        return CONTROL_STAGE_DESIGN_KEY_REQUIRED
    if "collision" in text or "duplicate high-confidence plasma-to-animal" in text:
        return PLASMA_ANIMAL_LINKAGE_COLLISION
    if "confidence provenance" in text or "confidence column" in text:
        return PLASMA_LINKAGE_CONFIDENCE_MISSING
    if "overlapping animals" in text or "linked animals" in text:
        return INSUFFICIENT_LINKED_ANIMALS
    if "animal_id" in text or "plasma-to-animal" in text:
        if "only" in text or "too few" in text or "minimum" in text or "overlap" in text:
            return INSUFFICIENT_LINKED_ANIMALS
        return PLASMA_BULK_ANIMAL_LINKAGE_MISSING
    if "plasma samples could be linked" in text or "high-confidence plasma" in text:
        return PLASMA_BULK_ANIMAL_LINKAGE_MISSING
    if "common tissues" in text or "shared tissue" in text or "overlapping tissues" in text:
        return INSUFFICIENT_SHARED_TISSUE_OVERLAP
    if "too few samples" in text or "insufficient samples" in text or "min_per_group" in text:
        return INSUFFICIENT_GROUP_SAMPLE_SIZE
    if (
        "unavailable" in text
        or "not available" in text
        or "missing effect tables" in text
        or "empty effect tables" in text
    ):
        return DATA_FILE_UNAVAILABLE
    if "not implemented" in text or "pseudobulk validation" in text:
        return SUBSET_PSEUDOBULK_NOT_IMPLEMENTED
    if "must contain" in text or "column missing" in text or "schema" in text:
        return SCHEMA_REQUIREMENT_MISSING
    if "not estimable" in text or "no sensitivity analyses" in text:
        return NO_ESTIMABLE_FEATURES
    if "failed" in text:
        return METHOD_RUNTIME_FAILURE
    if estimable and reason_l:
        return LOW_SAMPLE_SIZE_CAUTION
    return NON_ESTIMABLE_UNSPECIFIED


def missing_author_key_for_reason_code(reason_code: str) -> str:
    code = str(reason_code or "").strip() or NON_ESTIMABLE_UNSPECIFIED
    return MISSING_AUTHOR_KEY_BY_REASON_CODE.get(code, AUTHOR_KEY_REVIEW_REQUIRED)


def reason_fields_for_record(record: Mapping[str, Any]) -> dict[str, str]:
    existing_code = _clean_optional(record.get("reason_code"))
    existing_key = _clean_optional(record.get("missing_author_key"))
    estimable = _as_bool(record.get("estimable"), default=True)
    reason = _clean_optional(record.get("reason"))
    method = _clean_optional(record.get("method"))

    code = existing_code or infer_reason_code(
        reason=reason,
        method=method,
        estimable=estimable,
    )
    key = existing_key or missing_author_key_for_reason_code(code)
    if estimable and code == OK:
        key = ""
    return {"reason_code": code, "missing_author_key": key}


def annotate_reason_record(record: Mapping[str, Any]) -> dict[str, Any]:
    out = dict(record)
    out.update(reason_fields_for_record(out))
    return out


def annotate_reason_fields(
    df: pd.DataFrame,
    *,
    default_reason_code: Optional[str] = None,
    default_missing_author_key: Optional[str] = None,
) -> pd.DataFrame:
    out = df.copy()
    if "reason_code" not in out.columns:
        out["reason_code"] = ""
    if "missing_author_key" not in out.columns:
        out["missing_author_key"] = ""

    for idx, row in out.iterrows():
        fields = reason_fields_for_record(row.to_dict())
        if default_reason_code and not _clean_optional(out.at[idx, "reason_code"]):
            fields["reason_code"] = str(default_reason_code)
            fields["missing_author_key"] = missing_author_key_for_reason_code(
                fields["reason_code"]
            )
        if default_missing_author_key and not _clean_optional(out.at[idx, "missing_author_key"]):
            fields["missing_author_key"] = str(default_missing_author_key)
        out.at[idx, "reason_code"] = fields["reason_code"]
        out.at[idx, "missing_author_key"] = fields["missing_author_key"]

    return out
