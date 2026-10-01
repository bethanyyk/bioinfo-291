#!/usr/bin/env python
"""End-to-end driver for the GSE120575 pre/post-ICB analysis.

Each stage declares its output paths and is skipped when they all exist, so a
rerun after an interruption resumes rather than restarts. ``--force`` ignores
existing outputs; ``--stage`` runs one stage and its prerequisites' outputs are
assumed present.

    python workflows/run_pipeline.py --config configs/analysis.yaml
    python workflows/run_pipeline.py --stage abundance --force
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import anndata as ad  # noqa: E402

from icb_scrna import abundance as ab  # noqa: E402
from icb_scrna import annotate as an  # noqa: E402
from icb_scrna import de as DE  # noqa: E402
from icb_scrna import ingest  # noqa: E402
from icb_scrna import plots as PL  # noqa: E402
from icb_scrna import preprocess as pp  # noqa: E402
from icb_scrna import signature as sg  # noqa: E402
from icb_scrna.config import Config, load_config  # noqa: E402

CONTRASTS = {
    "treatment": ("timepoint", "Post", "Pre", None),
    "response_baseline": ("response", "Responder", "Non-responder", ("timepoint", "Pre")),
    "response_ontreatment": ("response", "Responder", "Non-responder", ("timepoint", "Post")),
}


# --------------------------------------------------------------------------
# stages
# --------------------------------------------------------------------------
def stage_download(cfg: Config) -> list[Path]:
    base = cfg["dataset"]["ftp_base"]
    manifest = []
    outs = []
    for fn in cfg["dataset"]["files"].values():
        dest = cfg.path("raw", fn)
        if not dest.exists():
            urllib.request.urlretrieve(base + fn, dest)
        h = hashlib.sha256()
        with open(dest, "rb") as fh:
            for blk in iter(lambda: fh.read(1 << 20), b""):
                h.update(blk)
        manifest.append({"file": fn, "url": base + fn,
                         "bytes": dest.stat().st_size, "sha256": h.hexdigest()})
        outs.append(dest)
    with open(cfg.path("metadata", "raw_download_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    return outs


def stage_ingest(cfg: Config) -> list[Path]:
    exp = cfg.path("raw", cfg["dataset"]["files"]["expression"])
    met = cfg.path("raw", cfg["dataset"]["files"]["cell_metadata"])
    meta = ingest.parse_cell_metadata(met)
    cells, header_labels = ingest.read_expression_headers(exp)
    if list(meta.index) != cells:
        raise ValueError("annotation row order does not match the expression header")
    meta = ingest.annotate_sort_protocol(meta, header_labels)
    X, genes = ingest.stream_expression(
        exp, n_cells=len(cells),
        chunk_rows=cfg["ingest"]["chunk_rows"],
        min_cells_per_gene=cfg["ingest"]["min_cells_per_gene"],
    )
    adata = ad.AnnData(X=X, obs=meta.copy(),
                       var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    adata.uns["accession"] = cfg["dataset"]["accession"]
    n_exp = cfg["dataset"]["expected_cells"]
    if adata.n_obs != n_exp:
        raise ValueError(f"expected {n_exp} cells, built {adata.n_obs}")
    meta.to_csv(cfg.path("metadata", "cell_annotation.csv"))
    ingest.sample_table(meta).to_csv(cfg.path("metadata", "sample_table.csv"), index=False)
    out = cfg.path("processed", "gse120575_raw.h5ad")
    adata.write_h5ad(out, compression="gzip")
    return [out]


def stage_atlas(cfg: Config) -> list[Path]:
    adata = ad.read_h5ad(cfg.path("processed", "gse120575_raw.h5ad"))
    adata = pp.add_qc_metrics(pp.to_natural_log(adata))
    qc_before = adata.obs[["n_genes", "pct_mito"]].copy()
    adata, tally = pp.filter_cells(adata, cfg["qc"])
    tally.to_csv(cfg.path("results", "qc_tally.csv"), index=False)
    adata = pp.build_atlas(adata, cfg["atlas"], seed=cfg.seed)
    pp.patient_entropy(adata).to_csv(cfg.path("results", "cluster_patient_entropy.csv"),
                                     index=False)
    sc.tl.leiden(adata, restrict_to=("leiden", ["0"]), resolution=0.4,
                 key_added="leiden_sub", random_state=cfg.seed,
                 flavor="igraph", n_iterations=2, directed=False)
    cols = an.score_panels(adata, seed=cfg.seed)
    adata = an.apply_labels(adata, an.CLUSTER_CELL_TYPE, "leiden_sub", "cell_type")
    adata = an.apply_labels(adata, an.CLUSTER_CELL_STATE, "leiden_sub", "cell_state")
    ev = an.cluster_panel_matrix(adata, cols, cluster_key="leiden_sub").round(3)
    ev.insert(0, "cell_type", pd.Series(an.CLUSTER_CELL_TYPE).reindex(ev.index).values)
    ev.to_csv(cfg.path("results", "cluster_annotation_evidence.csv"))
    fig = PL.qc_figure(qc_before, cfg["qc"])
    fig.savefig(cfg.path("figures", "fig0_qc.png"), dpi=300, bbox_inches="tight")
    out = cfg.path("processed", "gse120575_annotated.h5ad")
    adata.write_h5ad(out, compression="gzip")
    return [out]


def stage_abundance(cfg: Config) -> list[Path]:
    adata = _annotated(cfg)
    smeta = pd.read_csv(cfg.path("metadata", "sample_table.csv")).set_index("sample")
    frames = []
    for level in ("cell_type", "cell_state"):
        cnt, props = ab.proportion_table(
            adata.obs, cell_type_key=level,
            min_cells=cfg["abundance"]["min_cells_per_sample"])
        for name, (grp, test, ref, subset) in CONTRASTS.items():
            c = cnt
            if subset:
                c = c.loc[(smeta.loc[c.index, subset[0]].astype(str) == subset[1]).values]
            r = ab.test_abundance(c, smeta, grp, test, ref,
                                  covariates=["therapy"] if name == "treatment" else None,
                                  alpha=cfg["abundance"]["alpha"])
            r.insert(0, "level", level)
            r.insert(1, "analysis", name)
            frames.append(r)
        if level == "cell_type":
            ab.paired_change(props, smeta).to_csv(
                cfg.path("results", "differential_abundance_paired.csv"), index=False)
            PL.composition_figure(
                props, smeta,
                frames[-2],  # response_baseline for this level
                alpha=cfg["abundance"]["alpha"],
            ).savefig(cfg.path("figures", "fig2_composition.png"), dpi=300,
                      bbox_inches="tight")
    out = cfg.path("results", "differential_abundance.csv")
    pd.concat(frames, ignore_index=True).to_csv(out, index=False)
    return [out]


def stage_markers(cfg: Config) -> list[Path]:
    adata = _annotated(cfg)
    sc.tl.rank_genes_groups(adata, groupby="cell_type", method=cfg["markers"]["method"],
                            pts=True, key_added="markers_cell_type")
    M = sc.get.rank_genes_groups_df(adata, group=None, key="markers_cell_type").rename(
        columns={"group": "cell_type", "names": "gene", "logfoldchanges": "log2FC",
                 "pvals": "p_value", "pvals_adj": "p_adj_BH",
                 "pct_nz_group": "pct_in_group", "pct_nz_reference": "pct_in_rest"})
    full = cfg.path("results", "markers_cell_type_full.csv.gz")
    M.to_csv(full, index=False, compression="gzip")
    top = (M[(M["p_adj_BH"] < cfg["markers"]["max_padj"])
             & (M["log2FC"] >= cfg["markers"]["min_log2fc"])]
           .sort_values(["cell_type", "scores"], ascending=[True, False])
           .groupby("cell_type", observed=True).head(cfg["markers"]["n_top"]))
    out = cfg.path("results", "markers_cell_type_top.csv")
    top.to_csv(out, index=False)
    return [full, out]


def stage_de(cfg: Config) -> list[Path]:
    adata = _annotated(cfg)
    cts = list(adata.obs["cell_type"].cat.categories)
    pb = cfg["pseudobulk"]
    frames = []
    for name, (grp, test, ref, subset) in CONTRASTS.items():
        r = DE.run_contrast(adata, cts, grp, test, ref, subset,
                            min_cells=pb["min_cells"],
                            min_samples_per_group=pb["min_samples_per_group"])
        if not r.empty:
            r.insert(0, "analysis", name)
            frames.append(r)
    D = pd.concat(frames, ignore_index=True)
    full = cfg.path("results", "pseudobulk_de_full.csv.gz")
    D.to_csv(full, index=False, compression="gzip")
    sig = D[D["p_adj_BH"] < pb["max_padj"]]
    out = cfg.path("results", "pseudobulk_de_significant.csv")
    sig.to_csv(out, index=False)
    v = sig[(sig.analysis == "response_ontreatment") & (sig.cell_type == "CD8 T")]
    if not v.empty:
        sub = D[(D.analysis == "response_ontreatment") & (D.cell_type == "CD8 T")]
        PL.volcano_figure(sub, "CD8 T cells on treatment", "responders",
                          "non-responders", padj_cut=pb["max_padj"]).savefig(
            cfg.path("figures", "fig3_cd8_volcano.png"), dpi=300, bbox_inches="tight")
    return [full, out]


def stage_signature(cfg: Config) -> list[Path]:
    from icb_scrna.features import baseline_features, cd8_pseudobulk_by_patient
    adata = _annotated(cfg)
    smeta = pd.read_csv(cfg.path("metadata", "sample_table.csv")).set_index("sample")
    feats, y = baseline_features(adata, smeta)
    PB_pre, PB_post, y_post, pat_post = cd8_pseudobulk_by_patient(adata, smeta, y)
    s = cfg["signature"]
    results, scores = [], {}
    for name, F in feats.items():
        yy = y.loc[F.index]
        auc, sc_ = sg.loo_auc(F, yy, seed=cfg.seed)
        lo, hi = sg.bootstrap_auc_ci(yy, sc_, n_boot=s["n_bootstrap"], seed=cfg.seed)
        p = sg.permutation_p(auc, lambda ys, F=F: sg.loo_auc(F, ys, seed=cfg.seed)[0],
                             yy, n_perm=s["n_permutations"], seed=cfg.seed)
        results.append(dict(signature=name, n_patients=len(F), n_features=F.shape[1],
                            selection="fixed", auc_loo=auc, ci95_lo=lo, ci95_hi=hi,
                            permutation_p=p, genes_on_full_data=""))
        scores[name] = pd.Series(sc_, index=F.index)
    pat_pre = pd.Series(PB_pre.index, index=PB_pre.index)
    nested = {
        f"CD8 baseline gene program (nested, {s['n_genes']} genes)": {},
        f"CD8 on-treatment program scored at baseline (nested, {s['n_genes']} genes)":
            dict(selection_profiles=PB_post, selection_y=y_post, selection_patient=pat_post),
    }
    for name, kw in nested.items():
        yy = y.loc[PB_pre.index]
        auc, sc_, genes = sg.loo_auc_nested(PB_pre, yy, n_genes=s["n_genes"],
                                            patient=pat_pre, seed=cfg.seed, **kw)
        lo, hi = sg.bootstrap_auc_ci(yy, sc_, n_boot=s["n_bootstrap"], seed=cfg.seed)

        def run_cv(ys, kw=kw):
            k = dict(kw)
            if k.get("selection_profiles") is None:
                k["selection_y"] = ys
            return sg.loo_auc_nested(PB_pre, ys, n_genes=s["n_genes"], patient=pat_pre,
                                     seed=cfg.seed, **k)[0]
        p = sg.permutation_p(auc, run_cv, yy, n_perm=s["n_permutations"], seed=cfg.seed)
        results.append(dict(signature=name, n_patients=len(PB_pre), n_features=s["n_genes"],
                            selection="nested in CV", auc_loo=auc, ci95_lo=lo, ci95_hi=hi,
                            permutation_p=p, genes_on_full_data=";".join(genes)))
        scores[name] = pd.Series(sc_, index=PB_pre.index)
    perf = pd.DataFrame(results)
    out = cfg.path("results", "signature_performance.csv")
    perf.to_csv(out, index=False)
    S = pd.DataFrame(scores)
    S.insert(0, "responder", y.reindex(S.index).values)
    S.to_csv(cfg.path("results", "signature_scores_loo.csv"))
    PL.signature_figure(perf, scores, y).savefig(
        cfg.path("figures", "fig4_signatures.png"), dpi=300, bbox_inches="tight")
    return [out]


def stage_atlas_figure(cfg: Config) -> list[Path]:
    adata = _annotated(cfg)
    out = cfg.path("figures", "fig1_atlas.png")
    PL.atlas_figure(adata.obs, adata.obsm["X_umap"]).savefig(out, dpi=300,
                                                             bbox_inches="tight")
    return [out]


def _annotated(cfg: Config):
    return ad.read_h5ad(cfg.path("processed", "gse120575_annotated.h5ad"))


STAGES = {
    "download": (stage_download, ["data/raw/GSE120575_patient_ID_single_cells.txt.gz"]),
    "ingest": (stage_ingest, ["data/processed/gse120575_raw.h5ad"]),
    "atlas": (stage_atlas, ["data/processed/gse120575_annotated.h5ad"]),
    "atlas_figure": (stage_atlas_figure, ["figures/fig1_atlas.png"]),
    "abundance": (stage_abundance, ["results/differential_abundance.csv"]),
    "markers": (stage_markers, ["results/markers_cell_type_top.csv"]),
    "de": (stage_de, ["results/pseudobulk_de_significant.csv"]),
    "signature": (stage_signature, ["results/signature_performance.csv"]),
}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="configs/analysis.yaml")
    ap.add_argument("--stage", choices=list(STAGES), action="append")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    for name in (args.stage or list(STAGES)):
        fn, outputs = STAGES[name]
        paths = [cfg.root / o for o in outputs]
        if not args.force and all(p.exists() for p in paths):
            print(f"[skip] {name} - outputs present")
            continue
        print(f"[run ] {name}")
        produced = fn(cfg)
        print(f"[done] {name} -> {', '.join(str(p.relative_to(cfg.root)) for p in produced)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
