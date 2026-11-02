import shutil
from pathlib import Path

import multiqc
import pytest
from multiqc import report


def _run_module_on_dir(data_dir: Path):
    """Run the duplex_seq module with MultiQC file discovery on a temp directory."""
    multiqc.reset()
    multiqc.parse_logs(str(data_dir), run_modules=["duplex_seq"], quiet=True)
    assert len(report.modules) == 1
    return report.modules[0]


def test_parses_two_sample_csv_files(tmp_path):
    fixtures_dir = Path(__file__).parent
    shutil.copy(fixtures_dir / "NanoMB1Rep1_test_10k.csv", tmp_path / "NanoMB1Rep1_test_10k.csv")
    shutil.copy(fixtures_dir / "NanoMB2Rep1_test_10k.csv", tmp_path / "NanoMB2Rep1_test_10k.csv")

    module = _run_module_on_dir(tmp_path)

    assert "NanoMB1Rep1_test_10k" in module._sample_data
    assert "NanoMB2Rep1_test_10k" in module._sample_data
    assert module._sample_data["NanoMB1Rep1_test_10k"]["efficiency"] == pytest.approx(0.0547748866871889)
    assert module._sample_data["NanoMB2Rep1_test_10k"]["efficiency"] == pytest.approx(0.0512334521098765)


def test_parses_uppercase_headers(tmp_path):
    uppercase_csv = tmp_path / "uppercase_headers.csv"
    uppercase_csv.write_text(
        "\n".join(
            [
                "Sample,Metric,Value",
                "UpperSample,Efficiency,0.1",
                "UpperSample,drop_out_rate,0.2",
            ]
        )
    )

    module = _run_module_on_dir(tmp_path)

    assert "UpperSample" in module._sample_data
    assert module._sample_data["UpperSample"]["Efficiency"] == pytest.approx(0.1)
    assert module._sample_data["UpperSample"]["drop_out_rate"] == pytest.approx(0.2)


def test_grouped_section_anchors_are_unique(tmp_path):
    fixtures_dir = Path(__file__).parent
    shutil.copy(fixtures_dir / "NanoMB1Rep1_test_10k.csv", tmp_path / "NanoMB1Rep1_test_10k.csv")
    shutil.copy(fixtures_dir / "NanoMB2Rep1_test_10k.csv", tmp_path / "NanoMB2Rep1_test_10k.csv")

    module = _run_module_on_dir(tmp_path)

    anchors = [str(s.anchor) for s in module.sections]
    assert len(anchors) == len(set(anchors))
    assert "duplex_seq_on_target_rate" in anchors
    assert "duplex_seq_on_target_coverage" in anchors
    assert "duplex_seq_on_target_duplex_ratio_plot" in anchors
    assert "duplex_seq_gc_metrics" in anchors
    assert "duplex_seq_within_family_stats" in anchors
    assert "duplex_seq_family_size_comparison" in anchors
    assert anchors[-1] == "duplex_seq_frac_singletons_plot"


def test_groups_on_target_rate_raw_and_duplex(tmp_path):
    on_target_csv = tmp_path / "on_target_rates.csv"
    on_target_csv.write_text(
        "\n".join(
            [
                "sample,metric,value",
                "SampleA,efficiency,0.1",
                "SampleA,drop_out_rate,0.2",
                "SampleA,on_target_rate_raw,0.3",
                "SampleA,on_target_rate_duplex,0.4",
            ]
        )
    )

    module = _run_module_on_dir(tmp_path)

    anchors = [str(s.anchor) for s in module.sections]
    assert "duplex_seq_on_target_rate" in anchors
    assert "duplex_seq_on_target_rate_raw_plot" not in anchors
    assert "duplex_seq_on_target_rate_duplex_plot" not in anchors


def test_groups_on_target_coverage_and_plots_ratio(tmp_path):
    on_target_csv = tmp_path / "on_target_coverage.csv"
    on_target_csv.write_text(
        "\n".join(
            [
                "sample,metric,value",
                "SampleA,on_target_coverage_raw,176636.497683199",
                "SampleA,on_target_coverage_duplex,4147.76941198496",
                "SampleA,on_target_duplex_ratio,42.5859010322051",
            ]
        )
    )

    module = _run_module_on_dir(tmp_path)

    anchors = [str(s.anchor) for s in module.sections]
    assert "duplex_seq_on_target_coverage" in anchors
    assert "duplex_seq_on_target_coverage_raw_plot" not in anchors
    assert "duplex_seq_on_target_coverage_duplex_plot" not in anchors
    assert "duplex_seq_on_target_duplex_ratio_plot" in anchors
