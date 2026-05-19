import numpy as np
import pandas as pd
from pathlib import Path
import sys
import zipfile
from tempfile import TemporaryDirectory

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.exosome_effect import estimate_exosome_fraction_with_uncertainty
from src.attribution import build_mouse_exosome_metadata, compute_exosome_alignment_tables
from src.omix007582_audit import build_omix007582_sample_map_audit
from src.omix_io import load_omix_matrix
from src.repo_promotion import build_promotion_plan, apply_promotion_plan
from src.linkage_audit import build_estimability_report
from src.run_pipeline import _build_default_config, _causal_gate_reason, _ensure_standard_schema, STANDARD_RESULT_COLUMNS
from src.run_pipeline import _evidence_level


def _workspace_tempdir() -> TemporaryDirectory:
    base = Path(__file__).resolve().parents[1] / ".pytest_tmp"
    base.mkdir(exist_ok=True)
    return TemporaryDirectory(dir=base)


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


def test_mouse_exosome_parser_recovers_expected_fields():
    meta = build_mouse_exosome_metadata(["21_M_brain_GES_1", "8_W_liver_Ctrl_10"])
    first = meta.loc[meta["sample_id"] == "21_M_brain_GES_1"].iloc[0]
    second = meta.loc[meta["sample_id"] == "8_W_liver_Ctrl_10"].iloc[0]

    assert bool(first["parse_ok"]) is True
    assert float(first["age"]) == 21.0
    assert str(first["sex"]) == "M"
    assert str(first["tissue"]) == "brain"
    assert str(first["arm"]) == "GES"

    assert bool(second["parse_ok"]) is True
    assert str(second["sex"]) == "F"
    assert str(second["arm"]) == "Ctrl"


def test_exosome_alignment_unestimable_has_explicit_reason():
    prim = pd.DataFrame(
        {"mean_effect": [-0.4, -0.2]},
        index=["Hippocampus", "Liver_L"],
    )
    prim.index.name = "tissue"
    mouse = pd.DataFrame(
        [
            {
                "tissue": "brain",
                "contrast": "GES_vs_Veh",
                "mean_effect": -0.3,
                "estimable": True,
                "available": True,
                "ci_low": -0.5,
                "ci_high": -0.1,
                "permutation_p_value": 0.04,
            },
            {
                "tissue": "liver",
                "contrast": "WT_vs_Veh",
                "mean_effect": -0.1,
                "estimable": True,
                "available": True,
                "ci_low": -0.2,
                "ci_high": 0.0,
                "permutation_p_value": 0.2,
            },
        ]
    )
    by_tissue, summary = compute_exosome_alignment_tables(
        prim,
        mouse,
        tissue_map={"brain": "Hippocampus", "liver": "Liver_L"},
        contrasts=["GES_vs_Veh"],
        min_common_tissues=3,
        n_bootstrap=50,
        n_permutations=50,
        random_state=1,
    )
    assert bool(summary.loc[0, "estimable"]) is False
    assert isinstance(summary.loc[0, "reason"], str)
    assert len(summary.loc[0, "reason"]) > 0
    assert "common tissues" in summary.loc[0, "reason"]
    assert "reason" in by_tissue.columns


def test_evidence_level_ladder_supports_alignment_and_mediation():
    assert _evidence_level(estimable=False) == 0
    assert _evidence_level(estimable=True) == 1
    assert _evidence_level(estimable=True, has_linkage_support=True) == 2
    assert _evidence_level(estimable=True, has_exosome_alignment=True) == 3
    assert _evidence_level(
        estimable=True,
        tier="fully_linked",
        has_linked_mediation=True,
    ) == 4


def test_omix007582_audit_refuses_underdetermined_biological_mapping():
    with _workspace_tempdir() as tmpdir:
        tmp_path = Path(tmpdir)
        matrix_path = tmp_path / "OMIX007582_beta_matrix.csv"
        metadata_path = tmp_path / "OMIX007582-02.csv"
        idat_dir = tmp_path / "OMIX007582_idat"
        zip_path = tmp_path / "OMIX007582-03.zip"
        idat_dir.mkdir()

        pd.DataFrame(
            {
                "207925070004_R01C01": [0.1, 0.2],
                "207925070004_R01C02": [0.3, 0.4],
            },
            index=["cg1", "cg2"],
        ).to_csv(matrix_path)
        pd.DataFrame(
            {
                "OriginalSampleName": ["FV1-Heart", "FWT2-Heart"],
                "OriginalSampleName.1": ["FV1-Heart", "FWT2-Heart"],
                "sample": ["FV1", "FWT2"],
                "tissue": ["Heart", "Heart"],
                "group": ["O_V", "O_WT"],
            }
        ).to_csv(metadata_path, index=False)

        for suffix in ("Grn", "Red"):
            (idat_dir / f"207925070004_R01C01_{suffix}.idat").write_text("", encoding="utf-8")
            (idat_dir / f"207925070004_R01C02_{suffix}.idat").write_text("", encoding="utf-8")

        with zipfile.ZipFile(zip_path, "w") as handle:
            handle.writestr("207925070004_R01C01_Grn.idat", "")
            handle.writestr("207925070004_R01C01_Red.idat", "")
            handle.writestr("207925070004_R01C02_Grn.idat", "")
            handle.writestr("207925070004_R01C02_Red.idat", "")

        summary, overlap_audit, technical_inventory, metadata_inventory = build_omix007582_sample_map_audit(
            matrix_path=matrix_path,
            metadata_path=metadata_path,
            idat_dir=idat_dir,
            zip_path=zip_path,
        )

        assert bool(summary.loc[0, "available"]) is True
        assert bool(summary.loc[0, "estimable"]) is False
        assert str(summary.loc[0, "mapping_status"]) == "technical_ids_only"
        assert int(summary.loc[0, "exact_overlap_ids"]) == 0
        assert bool(summary.loc[0, "beta_equals_idat_directory"]) is True
        assert bool(summary.loc[0, "beta_equals_idat_archive"]) is True
        assert "cannot be recovered" in str(summary.loc[0, "reason"])
        assert int(overlap_audit["n_exact_overlaps"].sum()) == 0
        assert technical_inventory["mapped_sample_id"].isna().all()
        assert metadata_inventory["mapped_technical_id"].isna().all()


def test_build_default_config_supports_demo_profile():
    with _workspace_tempdir() as tmpdir:
        tmp_path = Path(tmpdir)
        processed = tmp_path / "data" / "PROCESSED"
        processed.mkdir(parents=True)
        (processed / "OMIX007580_01_example.txt").write_text("gene\ts1\nG1\t1\n", encoding="utf-8")
        pd.DataFrame({"sample_id": ["s1"], "group": ["Y_C"], "age": [4.0], "tissue": ["Brain"]}).to_csv(
            processed / "OMIX007580-02_example.csv",
            index=False,
        )
        pd.DataFrame({"Gene name": ["P1"], "FY_1": [1.0]}).to_csv(
            processed / "OMIX007581-01_example.csv",
            index=False,
        )
        pd.DataFrame({"s1": [0.1]}, index=["cg1"]).to_csv(processed / "OMIX007582_beta_matrix_example.csv")
        pd.DataFrame({"OriginalSampleName": ["s1"], "sample": ["s1"], "tissue": ["Brain"], "group": ["Y_C"]}).to_csv(
            processed / "OMIX007582-02_example.csv",
            index=False,
        )

        cfg = _build_default_config(tmp_path, profile="demo")
        assert cfg.data_profile == "demo"
        assert cfg.primate_bulk.matrix.name == "OMIX007580_01_example.txt"
        assert cfg.primate_plasma.matrix.name == "OMIX007581-01_example.csv"
        assert cfg.enable_mouse_exosome_block is False
        assert cfg.enable_subset_validation_block is False


def test_load_omix_matrix_accepts_string_feature_ids():
    with _workspace_tempdir() as tmpdir:
        tmp_path = Path(tmpdir)
        matrix_path = tmp_path / "OMIX007580_01_example.txt"
        matrix_path.write_text(
            "gene_id\ts1\ts2\nGene_0001\t1\t2\nGene_0002\t3\t4\n",
            encoding="utf-8",
        )

        matrix = load_omix_matrix(matrix_path, allowed_samples={"s1", "s2"})
        assert matrix.index.tolist() == ["Gene_0001", "Gene_0002"]
        assert matrix.columns.tolist() == ["s1", "s2"]
        assert float(matrix.loc["Gene_0001", "s1"]) == 1.0


def test_repo_promotion_plan_and_apply():
    with _workspace_tempdir() as tmpdir:
        tmp_path = Path(tmpdir)
        source_root = tmp_path / "source"
        target_root = tmp_path / "target"
        (source_root / "src").mkdir(parents=True)
        (target_root / "src").mkdir(parents=True)

        (source_root / "README.md").write_text("new readme\n", encoding="utf-8")
        (target_root / "README.md").write_text("old readme\n", encoding="utf-8")
        (source_root / "src" / "module.py").write_text("print('hello')\n", encoding="utf-8")

        manifest = {
            "sync_files": ["README.md", "src/module.py"],
            "optional_demo_files": [],
        }
        plan = build_promotion_plan(
            source_root=source_root,
            target_root=target_root,
            manifest=manifest,
            include_demo_data=False,
        )
        actions = {item.rel_path: item.action for item in plan}
        assert actions["README.md"] == "update"
        assert actions["src/module.py"] == "create"

        counts = apply_promotion_plan(plan)
        assert counts["update"] == 1
        assert counts["create"] == 1
        assert (target_root / "README.md").read_text(encoding="utf-8") == "new readme\n"
        assert (target_root / "src" / "module.py").read_text(encoding="utf-8") == "print('hello')\n"
