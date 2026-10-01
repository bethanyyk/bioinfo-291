"""Stage 1 - read the GEO supplementary files into an AnnData object.

The GSE120575 expression matrix is a 55k-gene x 16k-cell *dense tab-separated
text* file. Read naively it needs several GB of RAM. This module streams it in
row chunks, drops undetected genes per chunk, and accumulates sparse float32
blocks, so peak memory is set by ``ingest.chunk_rows`` rather than by the file.

File layout (verified against the downloaded files):

* expression: line 1 = cell barcodes (leading empty field), line 2 = sample
  labels (leading empty field), lines 3+ = ``gene<TAB>value * n_cells<TAB>``
  -- note the **trailing tab**, so data lines carry one more field than the
  header lines.
* annotation: 19 lines of GEO submission-template preamble, then a header row,
  then one row per cell.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

# Column header as spelled in the GEO file (the typo is theirs).
_SAMPLE_COL = "characteristics: patinet ID (Pre=baseline; Post= on treatment)"
_PREAMBLE_LINES = 19


def parse_cell_metadata(path: str | Path) -> pd.DataFrame:
    """Parse ``GSE120575_patient_ID_single_cells.txt.gz`` into a per-cell table.

    Returns a frame indexed by cell barcode with columns ``sample``, ``patient``,
    ``timepoint`` (``Pre``/``Post``), ``response``, ``therapy`` and
    ``enrichment`` (``CD45`` for the standard CD45+ sort, ``T_enriched`` for the
    T-cell-enriched biopsies, which were sorted differently and so carry a
    composition bias that must be modelled, not ignored).
    """
    raw = pd.read_csv(
        path,
        sep="\t",
        skiprows=_PREAMBLE_LINES,
        dtype=str,
        quoting=3,  # csv.QUOTE_NONE - the preamble contains unbalanced quotes
        encoding="latin-1",  # the protocol text carries a non-UTF-8 micro sign
    )
    raw = raw.dropna(subset=["title"])
    # the SAMPLES block ends at the first row that is not "Sample <n>"
    keep = raw["Sample name"].astype(str).str.match(r"Sample \d+$")
    raw = raw.loc[keep]

    meta = pd.DataFrame(index=pd.Index(raw["title"].values, name="cell"))
    meta["sample"] = raw[_SAMPLE_COL].values
    meta["response"] = raw["characteristics: response"].values
    meta["therapy"] = raw["characteristics: therapy"].values

    sample = meta["sample"].astype(str)
    meta["timepoint"] = np.where(sample.str.startswith("Post"), "Post", "Pre")
    # "Pre_P1", "Post_P1_2" -> patient P1
    meta["patient"] = sample.str.extract(r"^(?:Pre|Post)_(P\d+)", expand=False).values
    for col in ("sample", "response", "therapy", "timepoint", "patient"):
        meta[col] = meta[col].astype("category")
    return meta


def annotate_sort_protocol(meta: pd.DataFrame, header_labels: list[str]) -> pd.DataFrame:
    """Attach the sort protocol recorded in the expression file's second header row.

    The annotation file collapses every biopsy to ``Pre_Pn``/``Post_Pn``, but the
    expression header preserves suffixes such as ``_T_enriched`` and
    ``_myeloid_enriched``. Those cells were sorted on a lineage marker rather
    than on CD45 alone, so their cell-type proportions are set by the sort and
    carry no information about tumour composition. Composition analyses must
    restrict to ``enrichment == "CD45"``; expression analyses may keep them.
    """
    if len(header_labels) != len(meta):
        raise ValueError(
            f"{len(header_labels)} header labels for {len(meta)} cells"
        )
    labels = pd.Series(header_labels, index=meta.index, dtype=str)
    enrichment = np.where(
        labels.str.contains("T_enriched"),
        "T_enriched",
        np.where(labels.str.contains("myeloid_enriched"), "myeloid_enriched", "CD45"),
    )
    out = meta.copy()
    out["sort_protocol"] = labels.values
    out["enrichment"] = pd.Categorical(enrichment)
    return out


def sample_table(meta: pd.DataFrame) -> pd.DataFrame:
    """Collapse the per-cell table to one row per biopsy."""
    cols = ["patient", "timepoint", "response", "therapy", "enrichment"]
    tab = (
        meta.reset_index()
        .groupby("sample", observed=True)
        .agg(n_cells=("cell", "size"), **{c: (c, "first") for c in cols})
        .reset_index()
    )
    return tab.sort_values(["timepoint", "patient"]).reset_index(drop=True)


def read_expression_headers(path: str | Path) -> tuple[list[str], list[str]]:
    """Return (cell barcodes, per-cell sample label) from the two header lines."""
    with gzip.open(path, "rt") as fh:
        cells = fh.readline().rstrip("\n").split("\t")[1:]
        samples = fh.readline().rstrip("\n").split("\t")[1:]
    if len(cells) != len(samples):
        raise ValueError(
            f"header mismatch: {len(cells)} barcodes vs {len(samples)} sample labels"
        )
    return cells, samples


def stream_expression(
    path: str | Path,
    n_cells: int,
    chunk_rows: int = 2000,
    min_cells_per_gene: int = 10,
) -> tuple[sp.csr_matrix, list[str]]:
    """Stream the gene-by-cell matrix into a sparse cells-by-genes matrix.

    Genes detected in fewer than ``min_cells_per_gene`` cells are dropped inside
    the chunk loop, so they never contribute to the accumulated matrix.
    Returns ``(X, gene_names)`` with ``X`` of shape ``(n_cells, n_genes_kept)``.
    """
    blocks: list[sp.csr_matrix] = []
    genes: list[str] = []
    # Parse directly to float32. Letting pandas infer gives float64 and doubles
    # the per-chunk footprint, which matters on a 2 GB budget.
    dtype: dict[int, object] = {0: str}
    dtype.update({i: np.float32 for i in range(1, n_cells + 1)})
    reader = pd.read_csv(
        path,
        sep="\t",
        skiprows=2,
        header=None,
        index_col=0,
        usecols=range(n_cells + 1),  # drop the trailing empty field
        dtype=dtype,
        chunksize=chunk_rows,
        na_filter=False,
    )
    for chunk in reader:
        values = chunk.to_numpy(dtype=np.float32, copy=False)
        detected = (values > 0).sum(axis=1)
        keep = detected >= min_cells_per_gene
        if not keep.any():
            continue
        blocks.append(sp.csr_matrix(values[keep]))
        genes.extend(chunk.index[keep].astype(str).tolist())
        del values, chunk
    if not blocks:
        raise ValueError("no genes survived the detection filter")
    gene_by_cell = sp.vstack(blocks, format="csr")
    del blocks
    return gene_by_cell.T.tocsr(), genes
