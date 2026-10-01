"""Ingest tests. These build their own fixtures and never touch data/raw."""

from __future__ import annotations

import gzip

import numpy as np
import pandas as pd
import pytest

from icb_scrna import ingest

PREAMBLE = "\n".join(f"# preamble line {i}" for i in range(19))
HEADER = "\t".join(
    ["Sample name", "title", "source name", "organism",
     ingest._SAMPLE_COL, "characteristics: response", "characteristics: therapy"]
)


def _write_annotation(tmp_path, rows):
    p = tmp_path / "anno.txt.gz"
    body = "\n".join("\t".join(r) for r in rows)
    with gzip.open(p, "wt", encoding="latin-1") as fh:
        fh.write(f"{PREAMBLE}\n{HEADER}\n{body}\n")
    return p


def test_parse_cell_metadata_derives_patient_and_timepoint(tmp_path):
    rows = [
        ["Sample 1", "cellA", "Melanoma", "Homo sapiens", "Pre_P1", "Responder", "anti-PD1"],
        ["Sample 2", "cellB", "Melanoma", "Homo sapiens", "Post_P1_2", "Non-responder", "anti-PD1"],
        ["Sample 3", "cellC", "Melanoma", "Homo sapiens", "Post_P12", "Responder", "anti-CTLA4"],
    ]
    meta = ingest.parse_cell_metadata(_write_annotation(tmp_path, rows))
    assert list(meta.index) == ["cellA", "cellB", "cellC"]
    assert list(meta["timepoint"].astype(str)) == ["Pre", "Post", "Post"]
    assert list(meta["patient"].astype(str)) == ["P1", "P1", "P12"]


def test_parse_cell_metadata_stops_at_end_of_sample_block(tmp_path):
    rows = [
        ["Sample 1", "cellA", "Melanoma", "Homo sapiens", "Pre_P1", "Responder", "anti-PD1"],
        ["PROTOCOLS", "not a cell", "", "", "", "", ""],
        ["growth protocol", "tissue was disaggregated", "", "", "", "", ""],
    ]
    meta = ingest.parse_cell_metadata(_write_annotation(tmp_path, rows))
    assert list(meta.index) == ["cellA"]


def test_annotate_sort_protocol_flags_enriched_sorts():
    meta = pd.DataFrame(index=pd.Index(["a", "b", "c"], name="cell"))
    out = ingest.annotate_sort_protocol(
        meta, ["Post_P6_T_enriched", "Post_P17_myeloid_enriched", "Pre_P1"]
    )
    assert list(out["enrichment"].astype(str)) == ["T_enriched", "myeloid_enriched", "CD45"]


def test_annotate_sort_protocol_rejects_length_mismatch():
    meta = pd.DataFrame(index=["a", "b"])
    with pytest.raises(ValueError, match="header labels"):
        ingest.annotate_sort_protocol(meta, ["Pre_P1"])


def _write_expression(tmp_path, genes, values, n_cells):
    """Reproduce the GEO layout: two header rows, trailing tab on data rows."""
    p = tmp_path / "expr.txt.gz"
    cells = [f"cell{i}" for i in range(n_cells)]
    with gzip.open(p, "wt") as fh:
        fh.write("\t" + "\t".join(cells) + "\n")
        fh.write("\t" + "\t".join(["Pre_P1"] * n_cells) + "\n")
        for g, row in zip(genes, values):
            fh.write(g + "\t" + "\t".join(f"{v:.2f}" for v in row) + "\t\n")
    return p, cells


def test_read_expression_headers_strips_leading_field(tmp_path):
    p, cells = _write_expression(tmp_path, ["G1"], [[1.0, 2.0, 3.0]], 3)
    got_cells, got_samples = ingest.read_expression_headers(p)
    assert got_cells == cells
    assert got_samples == ["Pre_P1"] * 3


def test_stream_expression_is_cells_by_genes_and_filters_rare_genes(tmp_path):
    n_cells = 6
    genes = ["KEEP1", "RARE", "KEEP2"]
    values = [
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],   # detected in 6
        [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],   # detected in 1 -> dropped at min 3
        [0.0, 0.0, 7.0, 8.0, 9.0, 1.0],   # detected in 4
    ]
    p, _ = _write_expression(tmp_path, genes, values, n_cells)
    X, kept = ingest.stream_expression(p, n_cells=n_cells, chunk_rows=2,
                                       min_cells_per_gene=3)
    assert kept == ["KEEP1", "KEEP2"]
    assert X.shape == (n_cells, 2)
    np.testing.assert_allclose(X.toarray()[:, 0], values[0])
    np.testing.assert_allclose(X.toarray()[:, 1], values[2])


def test_stream_expression_chunk_size_does_not_change_result(tmp_path):
    rng = np.random.default_rng(0)
    n_cells, n_genes = 5, 11
    values = rng.random((n_genes, n_cells)) * (rng.random((n_genes, n_cells)) > 0.3)
    genes = [f"G{i}" for i in range(n_genes)]
    p, _ = _write_expression(tmp_path, genes, values, n_cells)
    a, ga = ingest.stream_expression(p, n_cells=n_cells, chunk_rows=2, min_cells_per_gene=1)
    b, gb = ingest.stream_expression(p, n_cells=n_cells, chunk_rows=100, min_cells_per_gene=1)
    assert ga == gb
    np.testing.assert_allclose(a.toarray(), b.toarray())
