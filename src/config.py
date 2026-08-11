"""
config.py

Typed configuration objects and defaults for the rejuvenation/exosome
pipeline, including dataset paths (OMIX IDs), column candidates, and
analysis hyperparameters used by run_pipeline.py.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional


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
    data_profile: str = "auto"
    data_root: Optional[Path] = None
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
    enable_methylation_block: bool = True
    enable_mouse_exosome_block: bool = True
    enable_subset_validation_block: bool = True
    enable_translation_module: bool = False

    # Mouse exosome alignment
    mouse_reference_arms: Optional[List[str]] = None
    mouse_contrasts: Optional[List[str]] = None
    mouse_min_reference_samples: int = 20
    mouse_min_reference_ages: int = 4
    mouse_min_samples_per_group: int = 3
    mouse_tissue_bootstrap: int = 1000
    mouse_tissue_permutations: int = 1000
    mouse_alignment_contrasts: Optional[List[str]] = None
    mouse_to_primate_tissue_map: Optional[Dict[str, str]] = None

    # Methylation validation
    methylation_group_age_map: Optional[Dict[str, float]] = None
    methylation_to_primate_tissue_map: Optional[Dict[str, str]] = None
    methylation_min_common_tissues: int = 3

    # Column candidates for auto-detection
    sample_id_col_candidates: Optional[List[str]] = None
    group_col_candidates: Optional[List[str]] = None
    tissue_col_candidates: Optional[List[str]] = None
    age_col_candidates: Optional[List[str]] = None
    sex_col_candidates: Optional[List[str]] = None
    animal_id_col_candidates: Optional[List[str]] = None

    # Labels
    mouse_treated_label: str = "GES"
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
            "Veh", "WT", "Ctrl", "Baseline"
        ]
        self.tissue_effect_covariates = self.tissue_effect_covariates or ["age", "sex", "batch"]
        self.sensitivity_top_feature_thresholds = self.sensitivity_top_feature_thresholds or [
            20,
            int(self.n_top_plasma_features),
            100,
        ]
        self.mouse_reference_arms = self.mouse_reference_arms or ["Baseline"]
        self.mouse_contrasts = self.mouse_contrasts or ["GES_vs_Veh", "WT_vs_Veh", "GES_vs_WT"]
        self.mouse_alignment_contrasts = self.mouse_alignment_contrasts or ["GES_vs_Veh", "WT_vs_Veh"]
        self.mouse_to_primate_tissue_map = self.mouse_to_primate_tissue_map or {
            "brain": "Hippocampus",
            "kidney": "Renal_cortex",
            "liver": "Liver_L",
            "lung": "Lung_R3_3",
            "muscle": "Quadriceps_muscle",
        }
        self.methylation_group_age_map = self.methylation_group_age_map or {
            "Y_C": 4.5,
            "Y_WT": 4.5,
            "Y_WS": 4.5,
            "M_C": 11.0,
            "O_C": 17.0,
            "O_V": 21.0,
            "O_WT": 21.0,
            "O_GES": 21.0,
            "O_CR": 21.0,
            "O_Met": 21.0,
            "O_VC": 21.0,
        }
        self.methylation_to_primate_tissue_map = self.methylation_to_primate_tissue_map or {
            "Brain": "Hippocampus",
            "Kidney": "Renal_cortex",
            "Liver": "Liver_L",
            "Lung": "Lung_R3_3",
            "Quadriceps muscle": "Quadriceps_muscle",
        }
