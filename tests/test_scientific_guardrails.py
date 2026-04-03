import numpy as np
import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.exosome_effect import estimate_exosome_fraction_with_uncertainty
from src.linkage_audit import build_estimability_report
from src.run_pipeline import _causal_gate_reason, _ensure_standard_schema, STANDARD_RESULT_COLUMNS


def test_unlinked_runs_do_not_pass_causal_gate():
    estimability = build_estimability_report(
        {
            "has_prim_animal_id_col": False,
            "has_plasma_animal_id_col": False,
            "n_overlap_animal_ids": 0,
        },
        min_samples_for_mediation=12,
    )
    reason = _causal_gate_reason(
        enable_mediation=True,
        enable_causal_decomposition=True,
        estimability_row=estimability,
    )
    assert reason is not None
    assert "animal_id" in reason.lower() or "overlap" in reason.lower()


def test_exosome_fraction_unestimable_has_explicit_reason():
    cells = pd.DataFrame({"mean_effect": [0.2, -0.1]}, index=["liver", "heart"])
    exo = pd.DataFrame(columns=["mean_effect"])
    out = estimate_exosome_fraction_with_uncertainty(
        effect_cells=cells,
        effect_exosomes=exo,
        min_common_tissues=3,
        n_bootstrap=100,
        n_permutations=100,
        random_state=1,
    )
    assert out["estimable"] is False
    assert isinstance(out["reason"], str)
    assert len(out["reason"]) > 0


def test_standard_schema_columns_always_present():
    df = pd.DataFrame([{"metric": 1.23}])
    out = _ensure_standard_schema(
        df,
        available=False,
        estimable=False,
        reason="stub",
        n_used=0,
        method="unit_test",
        ci_low=np.nan,
        ci_high=np.nan,
        evidence_level=0,
    )
    for col in STANDARD_RESULT_COLUMNS:
        assert col in out.columns
    assert out.loc[0, "reason"] == "stub"
