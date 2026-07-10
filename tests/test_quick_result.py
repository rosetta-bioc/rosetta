"""Tests for QuickResult class."""

import pytest
import pandas as pd
import numpy as np

from rosetta.quick_result import QuickResult
from rosetta.results import RosettaDataFrame


# --- Fixtures ---


@pytest.fixture
def seurat_result():
    """Simulated Seurat QuickResult."""
    clusters = pd.Series(
        [0, 0, 0, 1, 1, 1, 2, 2, 2, 2],
        index=[f"C{i}" for i in range(10)],
        name="seurat_clusters",
    )
    umap = pd.DataFrame(
        np.random.randn(10, 2),
        index=[f"C{i}" for i in range(10)],
        columns=["UMAP_1", "UMAP_2"],
    )
    data = {
        "clusters": clusters,
        "umap": umap,
        "variable_features": ["Gene1", "Gene2", "Gene3"],
    }
    metadata = {"pipeline": "standard", "n_pcs": 10}
    return QuickResult(data=data, method="seurat", metadata=metadata)


@pytest.fixture
def phyloseq_result():
    """Simulated phyloseq QuickResult."""
    diversity = pd.DataFrame(
        {"Shannon": [2.1, 2.5, 1.8, 2.3], "Simpson": [0.8, 0.85, 0.75, 0.82]},
        index=["S1", "S2", "S3", "S4"],
    )
    data = {"diversity": diversity, "measures": ["Shannon", "Simpson"]}
    metadata = {"measures": ["Shannon", "Simpson"]}
    return QuickResult(data=data, method="phyloseq", metadata=metadata)


@pytest.fixture
def deseq2_quick_result():
    """Simulated DESeq2 QuickResult."""
    results_df = RosettaDataFrame(
        {
            "baseMean": [100.0, 200.0, 300.0],
            "log2FoldChange": [1.5, -2.0, 0.5],
            "lfcSE": [0.3, 0.4, 0.5],
            "pvalue": [0.001, 0.01, 0.2],
            "padj": [0.01, 0.05, 0.4],
        },
        index=["GeneA", "GeneB", "GeneC"],
    )
    data = {"results": results_df}
    metadata = {"design": "~ condition", "alpha": 0.05}
    return QuickResult(data=data, method="deseq2", metadata=metadata)


@pytest.fixture
def edger_quick_result():
    """Simulated edgeR QuickResult."""
    results_df = RosettaDataFrame(
        {
            "logFC": [1.2, -1.8, 0.3],
            "logCPM": [5.0, 6.0, 4.5],
            "F": [12.0, 15.0, 2.0],
            "PValue": [0.001, 0.005, 0.3],
            "FDR": [0.01, 0.03, 0.5],
        },
        index=["GeneA", "GeneB", "GeneC"],
    )
    data = {"results": results_df}
    metadata = {"design": "~ condition"}
    return QuickResult(data=data, method="edger", metadata=metadata)


# --- Dict-like access tests ---


class TestDictAccess:
    def test_getitem(self, seurat_result):
        """Test that result['key'] works."""
        clusters = seurat_result["clusters"]
        assert isinstance(clusters, pd.Series)
        assert len(clusters) == 10

    def test_keys(self, seurat_result):
        """Test that result.keys() returns expected keys."""
        keys = set(seurat_result.keys())
        assert keys == {"clusters", "umap", "variable_features"}

    def test_contains(self, seurat_result):
        """Test that 'key' in result works."""
        assert "clusters" in seurat_result
        assert "umap" in seurat_result
        assert "nonexistent" not in seurat_result

    def test_len(self, seurat_result):
        """Test that len(result) returns number of data entries."""
        assert len(seurat_result) == 3

    def test_iter(self, seurat_result):
        """Test that iteration over result works."""
        keys = list(seurat_result)
        assert set(keys) == {"clusters", "umap", "variable_features"}

    def test_get_with_default(self, seurat_result):
        """Test get() with default value."""
        assert seurat_result.get("nonexistent", "default") == "default"
        assert seurat_result.get("clusters") is not None

    def test_values(self, seurat_result):
        """Test that result.values() works."""
        vals = list(seurat_result.values())
        assert len(vals) == 3

    def test_items(self, seurat_result):
        """Test that result.items() works."""
        items = dict(seurat_result.items())
        assert "clusters" in items

    def test_keyerror_on_missing(self, seurat_result):
        """Test that missing key raises KeyError."""
        with pytest.raises(KeyError):
            _ = seurat_result["nonexistent"]


# --- Property tests ---


class TestProperties:
    def test_data_property(self, seurat_result):
        assert isinstance(seurat_result.data, dict)

    def test_method_property(self, seurat_result):
        assert seurat_result.method == "seurat"

    def test_metadata_property(self, seurat_result):
        assert seurat_result.metadata == {"pipeline": "standard", "n_pcs": 10}

    def test_repr(self, seurat_result):
        r = repr(seurat_result)
        assert "QuickResult" in r
        assert "seurat" in r


# --- Report tests ---


class TestReport:
    def test_report_returns_string(self, seurat_result):
        """Test that .report() returns a string."""
        text = seurat_result.report()
        assert isinstance(text, str)
        assert len(text) > 0

    def test_report_prints_output(self, seurat_result, capsys):
        """Test that .report() prints to stdout."""
        seurat_result.report()
        captured = capsys.readouterr()
        assert len(captured.out) > 0


class TestSeuratReport:
    def test_includes_cluster_info(self, seurat_result):
        """Test that Seurat report includes cluster information."""
        text = seurat_result.report()
        assert "Seurat" in text
        assert "Cluster" in text or "cluster" in text
        assert "10" in text  # total cells
        assert "3" in text  # 3 clusters

    def test_includes_cell_count(self, seurat_result):
        """Test that Seurat report includes cell count."""
        text = seurat_result.report()
        assert "Total cells" in text

    def test_includes_cluster_sizes(self, seurat_result):
        """Test that Seurat report includes per-cluster sizes."""
        text = seurat_result.report()
        # Cluster 0: 3 cells, Cluster 1: 3 cells, Cluster 2: 4 cells
        assert "3" in text
        assert "4" in text


class TestPhyloseqReport:
    def test_includes_diversity_info(self, phyloseq_result):
        """Test that phyloseq report includes diversity metrics."""
        text = phyloseq_result.report()
        assert "Diversity" in text or "diversity" in text
        assert "Shannon" in text

    def test_includes_sample_count(self, phyloseq_result):
        """Test that phyloseq report includes sample count."""
        text = phyloseq_result.report()
        assert "4" in text  # 4 samples

    def test_includes_statistics(self, phyloseq_result):
        """Test that phyloseq report includes mean/sd/range."""
        text = phyloseq_result.report()
        assert "mean" in text
        assert "sd" in text
        assert "range" in text


class TestDeseq2Report:
    def test_delegates_to_rosettadataframe(self, deseq2_quick_result):
        """Test that DESeq2 report delegates to RosettaDataFrame.report."""
        text = deseq2_quick_result.report()
        assert "DESeq2" in text
        assert "Significant" in text

    def test_includes_gene_counts(self, deseq2_quick_result):
        """Test DESeq2 report includes gene counts."""
        text = deseq2_quick_result.report()
        assert "3" in text  # total genes


class TestEdgerReport:
    def test_delegates_to_rosettadataframe(self, edger_quick_result):
        """Test that edgeR report delegates to RosettaDataFrame.report."""
        text = edger_quick_result.report()
        assert "edgeR" in text
        assert "Significant" in text

    def test_includes_gene_counts(self, edger_quick_result):
        """Test edgeR report includes gene counts."""
        text = edger_quick_result.report()
        assert "3" in text  # total genes


class TestGenericReport:
    def test_unknown_method_uses_generic(self):
        """Test that an unknown method uses the generic reporter."""
        result = QuickResult(data={"key": "value"}, method="unknown")
        text = result.report()
        assert "QuickResult" in text
        assert "key" in text


class TestEmptyMetadata:
    def test_default_metadata_is_empty_dict(self):
        """Test that metadata defaults to empty dict."""
        result = QuickResult(data={"x": 1}, method="test")
        assert result.metadata == {}
