"""
config.py

Typed configuration objects and defaults for the rejuvenation/exosome
pipeline, including dataset paths (OMIX IDs), column candidates, and
analysis hyperparameters used by run_pipeline.py.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional


@dataclass(frozen=True)
class OmixPaths:
    """Container for OMIX matrix + metadata paths."""
    matrix: Path
    metadata: Optional[Path] = None


@dataclass
class PipelineConfig:
    """
    Configuration for the lightweight SRSC rejuvenation/exosome pipeline.

    This pipeline is designed for processed OMIX matrices (counts/beta/proteins)
    and avoids raw FASTQ/BAM downloads.

    The goal is mechanistic and translational inference, not full raw reproducibility.
    """

    # Core OMIX
    primate_bulk: OmixPaths
    primate_bulk_age_col = "age"
    primate_plasma: Optional[OmixPaths] = None
    primate_methylation: Optional[OmixPaths] = None
    mouse_exosome_bulk: Optional[OmixPaths] = None
    max_allowed_samples: int = 5000
    min_samples_for_mediation: int = 12
    min_samples_per_group_for_rejuv: int = 2
    mediation_bootstrap: int = 500
    clock_model: str = "ridge"  # or whatever you implemented
    control_label: str = "Control"
    primate_treated_label: str = "O_GES"
    primate_control_labels: Optional[List[str]] = None
    top_genes_per_tissue: int = 100
    top_plasma_biomarkers: int = 50
    random_seed: int = 42
    focus_genes: list[str] = ("FOXO3", "SRC")
    tissue_weighting: str = "uniform"
    exosome_fraction_method: str = "correlation_ratio"
    exosome_fraction_bootstrap: int = 2000
    exosome_fraction_permutations: int = 1000
    exosome_min_common_tissues: int = 3
    exosome_min_cells_median_abs: float = 1e-8
    enable_mediation: bool = True
    enable_causal_decomposition: bool = True
    # Optional modules flags
    enable_methylation_block: bool = False
    enable_translation_module: bool = False

    # Column candidates for auto-detection
    sample_id_col_candidates: Optional[List[str]] = None
    group_col_candidates: Optional[List[str]] = None
    tissue_col_candidates: Optional[List[str]] = None
    age_col_candidates: Optional[List[str]] = None
    sex_col_candidates: Optional[List[str]] = None
    animal_id_col_candidates: Optional[List[str]] = None

    # Labels
    mouse_treated_label: str = "Exosome"  # adjust to your metadata
    mouse_control_labels: Optional[List[str]] = None

    # Feature selection
    n_top_features_expr: int = 0
    n_top_features_clock: int = 2000
    n_top_features_proxy: int = 1000

    # Plasma score
    n_top_plasma_features: int = 50
    sensitivity_top_feature_thresholds: Optional[List[int]] = None
    plasma_biomarker_min_pairs: int = 8
    plasma_biomarker_bootstrap: int = 120
    plasma_biomarker_sign_agreement_min: float = 0.8
    plasma_biomarker_stability_top_k: int = 250

    # Statistics
    random_state: int = 42
    n_bootstrap: int = 2000

    # Output
    results_dir: Path = Path("../results")
    figures_dir: Path = Path("../figures")

    # Optional gene sets
    gmt_path: Optional[Path] = None

    # Plasma ranking plot
    plasma_ranking_requires_spearman: bool = False
    plasma_ranking_fallback_metric: str = "variance"  # "variance" | "abs_mean" | "abs_diff"
    tissue_effect_covariates: Optional[List[str]] = None

    def __post_init__(self):
        self.sample_id_col_candidates = self.sample_id_col_candidates or [
            "sample_id", "Sample", "sample", "SampleID", "sample_name"
        ]
        self.group_col_candidates = self.group_col_candidates or [
            "group", "Group", "treatment", "Treatment", "condition", "Condition"
        ]
        self.tissue_col_candidates = self.tissue_col_candidates or [
            "tissue", "Tissue", "organ", "Organ"
        ]
        self.age_col_candidates = self.age_col_candidates or [
            "agenumb",
            "age_num",
            "age_years",
            "Age (years)",
            "Age(years)",
            "age",
            "Age",
            "chrono_age",
            "chronological_age",
        ]
        self.sex_col_candidates = self.sex_col_candidates or [
            "sex", "Sex", "gender", "Gender"
        ]
        self.animal_id_col_candidates = self.animal_id_col_candidates or [
            "animal_id",
            "AnimalID",
            "donor_id",
            "Donor",
            "subject_id",
            "Subject",
            "orig.ident",
            "orig_ident",
            "orig.ident.id",
        ]

        self.primate_control_labels = self.primate_control_labels or [
            "Y_C", "M_C", "O_C", "O_WT", "O_V"
        ]
        self.mouse_control_labels = self.mouse_control_labels or [
            "Saline", "Control", "WTC"
        ]
        self.tissue_effect_covariates = self.tissue_effect_covariates or ["age", "sex", "batch"]
        self.sensitivity_top_feature_thresholds = self.sensitivity_top_feature_thresholds or [
            20,
            int(self.n_top_plasma_features),
            100,
        ]
