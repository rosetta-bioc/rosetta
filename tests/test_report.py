"""Tests for RosettaDataFrame.report() method."""

import numpy as np
import pandas as pd
import pytest

from rosetta.results import RosettaDataFrame


@pytest.fixture
def deseq2_results():
    np.random.seed(42)
    n = 1000
    return RosettaDataFrame({
        "baseMean": np.random.exponential(100, n),
        "log2FoldChange": np.random.normal(0, 1.5, n),
        "lfcSE": np.abs(np.random.normal(0.3, 0.1, n)),
        "stat": np.random.normal(0, 2, n),
        "pvalue": np.random.uniform(0, 1, n),
        "padj": np.concatenate([np.random.uniform(0, 0.01, 80), np.random.uniform(0.1, 1, n - 80)]),
    })


@pytest.fixture
def edger_results():
    np.random.seed(42)
    n = 500
    return RosettaDataFrame({
        "logFC": np.random.normal(0, 2, n),
        "logCPM": np.random.normal(5, 1, n),
        "F": np.abs(np.random.normal(3, 2, n)),
        "PValue": np.random.uniform(0, 1, n),
        "FDR": np.concatenate([np.random.uniform(0, 0.02, 40), np.random.uniform(0.1, 1, n - 40)]),
    })


@pytest.fixture
def enrichment_results():
    return RosettaDataFrame({
        "Description": ["cell cycle", "apoptosis", "DNA repair", "metabolism", "signaling"],
        "GeneRatio": ["20/200", "15/200", "10/200", "8/200", "5/200"],
        "p.adjust": [0.001, 0.01, 0.03, 0.08, 0.2],
        "Count": [20, 15, 10, 8, 5],
    })


class TestReport:
    def test_deseq2_report(self, deseq2_results):
        text = deseq2_results.report()
        assert "DESeq2 Results Summary" in text
        assert "1,000" in text
        assert "Upregulated" in text
        assert "Downregulated" in text

    def test_edger_report(self, edger_results):
        text = edger_results.report()
        assert "edgeR Results Summary" in text
        assert "500" in text

    def test_enrichment_report(self, enrichment_results):
        text = enrichment_results.report()
        assert "Enrichment Results Summary" in text
        assert "cell cycle" in text

    def test_custom_alpha(self, deseq2_results):
        strict = deseq2_results.report(alpha=0.001)
        lenient = deseq2_results.report(alpha=0.1)
        # Stricter alpha should show fewer significant genes
        assert "Significant" in strict
        assert "Significant" in lenient

    def test_dataframe_operations_preserved(self, deseq2_results):
        # Slicing should still work
        subset = deseq2_results[deseq2_results["padj"] < 0.05]
        assert isinstance(subset, RosettaDataFrame)
        assert hasattr(subset, "report")

    def test_to_csv_works(self, deseq2_results, tmp_path):
        path = tmp_path / "results.csv"
        deseq2_results.to_csv(path)
        reloaded = pd.read_csv(path, index_col=0)
        assert len(reloaded) == len(deseq2_results)


class TestMethodMetadata:
    def test_report_shows_method_when_set(self, deseq2_results):
        deseq2_results._rosetta_method = "deseq2"
        text = deseq2_results.report()
        assert "Method: deseq2" in text

    def test_report_no_method_line_when_unset(self, deseq2_results):
        # _rosetta_method not set — no "Method:" prefix
        text = deseq2_results.report()
        assert "Method:" not in text

    def test_method_line_appears_before_summary(self, deseq2_results):
        deseq2_results._rosetta_method = "deseq2"
        text = deseq2_results.report()
        method_pos = text.index("Method:")
        summary_pos = text.index("DESeq2 Results Summary")
        assert method_pos < summary_pos

    def test_method_propagates_through_slice(self, deseq2_results):
        deseq2_results._rosetta_method = "deseq2"
        subset = deseq2_results[deseq2_results["padj"] < 0.05]
        assert subset._rosetta_method == "deseq2"

    def test_edger_method_label(self, edger_results):
        edger_results._rosetta_method = "edger"
        text = edger_results.report()
        assert "Method: edger" in text


class TestPlotMethods:
    """Tests for .volcano() and .ma_plot() on RosettaDataFrame."""

    def test_volcano_returns_figure(self, deseq2_results):
        import matplotlib
        matplotlib.use("Agg")
        fig = deseq2_results.volcano()
        import matplotlib.pyplot as plt
        assert hasattr(fig, "savefig")
        plt.close("all")

    def test_volcano_with_highlight(self, deseq2_results):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        genes = list(deseq2_results.index[:2])
        fig = deseq2_results.volcano(highlight_genes=genes)
        assert fig is not None
        plt.close("all")

    def test_volcano_custom_alpha_and_lfc(self, deseq2_results):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig = deseq2_results.volcano(alpha=0.01, lfc_cutoff=2.0)
        assert fig is not None
        plt.close("all")

    def test_ma_plot_returns_figure(self, deseq2_results):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig = deseq2_results.ma_plot()
        assert hasattr(fig, "savefig")
        plt.close("all")

    def test_ma_plot_custom_title(self, deseq2_results):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig = deseq2_results.ma_plot(title="My MA Plot")
        ax = fig.axes[0]
        assert ax.get_title() == "My MA Plot"
        plt.close("all")

    def test_ma_plot_no_mean_col_raises(self):
        import matplotlib
        matplotlib.use("Agg")
        df = RosettaDataFrame({
            "log2FoldChange": [1.0, -1.0],
            "padj": [0.01, 0.5],
        })
        with pytest.raises(ValueError, match="mean expression"):
            df.ma_plot()

    def test_volcano_edger(self, edger_results):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig = edger_results.volcano()
        assert fig is not None
        plt.close("all")
