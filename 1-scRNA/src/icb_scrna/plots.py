"""Figure builders. Each returns a Matplotlib figure; the caller saves it.

Conventions enforced here: one colour per cell type threaded across every
panel, embeddings drawn without ticks, distributions shown as raw points when
n is small, and no claim in a title that the plotted data does not support.
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

# One hue per cell type, reused in every panel of every figure.
CELL_TYPE_COLORS: dict[str, str] = {
    "CD8 T": "#1f4e79", "CD4 T": "#4a90c2", "Treg": "#7fb3d5", "gdT": "#2e7d6f",
    "Cycling T": "#9467bd", "NK": "#17becf", "B": "#d95f02", "Plasma": "#8c4a0a",
    "Monocyte/Macrophage": "#b8860b", "cDC": "#e6ab02", "pDC": "#a6761d",
}
RESPONSE_COLORS = {"Responder": "#0072B2", "Non-responder": "#D55E00"}
TIMEPOINT_COLORS = {"Pre": "#7f7f7f", "Post": "#111111"}


def _embedding_axes(ax, xy, label="UMAP"):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    x0, y0 = xy[:, 0].min(), xy[:, 1].min()
    dx = 0.18 * (xy[:, 0].max() - x0); dy = 0.18 * (xy[:, 1].max() - y0)
    ax.annotate("", xy=(x0 + dx, y0 - dy * 0.35), xytext=(x0, y0 - dy * 0.35),
                arrowprops=dict(arrowstyle="->", lw=0.8, color="#444444"))
    ax.annotate("", xy=(x0 - dx * 0.35, y0 + dy), xytext=(x0 - dx * 0.35, y0),
                arrowprops=dict(arrowstyle="->", lw=0.8, color="#444444"))
    ax.text(x0 + dx * 0.5, y0 - dy * 0.75, f"{label} 1", ha="center", va="top", fontsize=6,
            color="#444444")
    ax.text(x0 - dx * 0.6, y0 + dy * 0.5, f"{label} 2", ha="right", va="center", fontsize=6,
            color="#444444", rotation=90)
    ax.margins(0.1)


def atlas_figure(obs: pd.DataFrame, umap: np.ndarray) -> plt.Figure:
    """Three embeddings: identity, timepoint, response."""
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.9))
    order = [c for c in CELL_TYPE_COLORS if c in set(obs["cell_type"].astype(str))]

    ax = axes[0]
    for ct in order:
        m = (obs["cell_type"].astype(str) == ct).to_numpy()
        ax.scatter(umap[m, 0], umap[m, 1], s=1.6, lw=0, alpha=0.75,
                   c=CELL_TYPE_COLORS[ct], rasterized=True)
    for ct in order:
        m = (obs["cell_type"].astype(str) == ct).to_numpy()
        cx, cy = np.median(umap[m, 0]), np.median(umap[m, 1])
        ax.text(cx, cy, ct, fontsize=6.5, ha="center", va="center", zorder=5,
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.82))
    ax.set_title("Eleven immune populations in CD45+ biopsies")
    _embedding_axes(ax, umap)

    for ax, key, colors, title in [
        (axes[1], "timepoint", TIMEPOINT_COLORS, "Timepoints overlap - no global shift on treatment"),
        (axes[2], "response", RESPONSE_COLORS, "Responders and non-responders overlap"),
    ]:
        for lev, col in colors.items():
            m = (obs[key].astype(str) == lev).to_numpy()
            ax.scatter(umap[m, 0], umap[m, 1], s=1.6, lw=0, alpha=0.6, c=col,
                       rasterized=True, label=f"{lev} (n={int(m.sum()):,} cells)")
        ax.set_title(title)
        _embedding_axes(ax, umap)
        leg = ax.legend(frameon=False, markerscale=6, fontsize=6.5, loc="lower right",
                        handletextpad=0.3, borderpad=0.1)
        for h in leg.legend_handles:
            h.set_alpha(1.0)
    fig.tight_layout()
    return fig


def composition_figure(props: pd.DataFrame, meta: pd.DataFrame,
                       stats_baseline: pd.DataFrame, alpha: float = 0.05) -> plt.Figure:
    """Per-biopsy composition and the baseline responder contrast."""
    fig = plt.figure(figsize=(12.5, 5.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.3, 1.0], wspace=0.30,
                          bottom=0.22, top=0.90)

    # -- stacked composition, biopsies ordered by group
    ax = fig.add_subplot(gs[0, 0])
    m = meta.loc[props.index]
    order = m.sort_values(["timepoint", "response"]).index
    P = props.loc[order]
    cts = [c for c in CELL_TYPE_COLORS if c in P.columns]
    bottom = np.zeros(len(P))
    for ct in cts:
        ax.bar(np.arange(len(P)), P[ct].to_numpy(), bottom=bottom, width=0.92,
               color=CELL_TYPE_COLORS[ct], lw=0, label=ct)
        bottom += P[ct].to_numpy()
    ax.set_xlim(-0.6, len(P) - 0.4); ax.set_ylim(0, 1)
    ax.set_ylabel("fraction of CD45+ cells in biopsy")
    ax.set_xticks([])
    groups = m.loc[order].apply(lambda r: f"{r['timepoint']}\n{r['response'][:3]}", axis=1)
    start = 0
    for lab in groups.unique():
        n = int((groups == lab).sum())
        ax.annotate(lab.replace("Non", "NR").replace("Res", "R"),
                    xy=(start + n / 2 - 0.5, -0.055), xycoords=("data", "axes fraction"),
                    ha="center", va="top", fontsize=6.5)
        if start > 0:
            ax.axvline(start - 0.5, color="white", lw=1.4)
        start += n
    ax.set_title("Composition of each biopsy")
    ax.legend(frameon=False, fontsize=6.2, ncol=6, loc="upper center",
              bbox_to_anchor=(0.5, -0.11), handlelength=1.0, handletextpad=0.4,
              columnspacing=1.0)

    # -- baseline responder vs non-responder, raw biopsies + FDR marks
    ax2 = fig.add_subplot(gs[0, 1])
    base = m[m["timepoint"].astype(str) == "Pre"].index
    st = stats_baseline.set_index("cell_type")
    rank = st["p_adj_BH"].reindex(cts).sort_values()
    ypos = np.arange(len(rank))[::-1]
    for yi, ct in zip(ypos, rank.index):
        for resp, dx in [("Responder", 0.17), ("Non-responder", -0.17)]:
            v = props.loc[[i for i in base if m.loc[i, "response"] == resp], ct].to_numpy()
            ax2.scatter(v, np.full(len(v), yi + dx), s=11, alpha=0.75, lw=0,
                        color=RESPONSE_COLORS[resp])
            ax2.plot([np.median(v)] * 2, [yi + dx - 0.1, yi + dx + 0.1],
                     color=RESPONSE_COLORS[resp], lw=1.6)
        q = st.loc[ct, "p_adj_BH"]
        if q < alpha:
            ax2.text(0.72, yi, f"q={q:.3f}", transform=ax2.get_yaxis_transform(),
                     fontsize=6, va="center", color="#222222")
    ax2.set_yticks(ypos); ax2.set_yticklabels(rank.index, fontsize=7)
    ax2.set_xlabel("fraction of CD45+ cells (baseline biopsies)")
    ax2.set_xlim(-0.02, 0.9)
    ax2.set_title("Higher B and CD4 T in responders at baseline")
    handles = [Line2D([], [], marker="o", ls="", color=c, label=k, markersize=4)
               for k, c in RESPONSE_COLORS.items()]
    ax2.legend(handles=handles, frameon=False, fontsize=6.5, loc="lower right")
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    return fig


def volcano_figure(de: pd.DataFrame, title: str, up_label: str, down_label: str,
                   padj_cut: float = 0.1, n_label: int = 14) -> plt.Figure:
    """Volcano for one cell type x contrast."""
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    fig.subplots_adjust(left=0.34, right=0.97)
    x = de["log2FC"].to_numpy()
    yv = -np.log10(de["p_value"].to_numpy())
    sig = (de["p_adj_BH"] < padj_cut).to_numpy()
    ax.scatter(x[~sig], yv[~sig], s=5, lw=0, color="#cccccc", rasterized=True)
    ax.scatter(x[sig & (x > 0)], yv[sig & (x > 0)], s=9, lw=0,
               color=RESPONSE_COLORS["Responder"])
    ax.scatter(x[sig & (x < 0)], yv[sig & (x < 0)], s=9, lw=0,
               color=RESPONSE_COLORS["Non-responder"])
    # Label the strongest hits in a spread column on each side rather than
    # beside each point: at this density inline labels collide.
    lab = de[sig].reindex(de[sig]["p_value"].sort_values().index).head(n_label)
    ax.margins(0.10)
    for side, sel in (("left", lab[lab["log2FC"] < 0]), ("right", lab[lab["log2FC"] > 0])):
        sel = sel.sort_values("p_value")
        if sel.empty:
            continue
        n = len(sel)
        # axes-fraction placement: evenly spaced by construction, so the
        # columns cannot collide however dense the point cloud is
        ys = np.linspace(0.97, 0.03, n) if n > 1 else [0.5]
        tx, ha = (-0.08, "right") if side == "left" else (1.03, "left")
        for (r, ty) in zip(sel.itertuples(), ys):
            ax.annotate(r.gene, xy=(r.log2FC, -np.log10(r.p_value)),
                        xytext=(tx, ty), textcoords="axes fraction",
                        fontsize=6, ha=ha, va="center", style="italic",
                        annotation_clip=False,
                        arrowprops=dict(arrowstyle="-", lw=0.35, color="#bbbbbb",
                                        shrinkA=1, shrinkB=1))
    ax.set_xlabel(f"log2 fold change   ({down_label} <-- --> {up_label})")
    ax.set_ylabel("-log10 p (Welch, pseudobulk)")
    ax.yaxis.set_label_coords(-0.30, 0.5)  # clear of the gene-label column
    ax.set_title(title)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return fig


def signature_figure(perf: pd.DataFrame, scores: dict[str, pd.Series],
                     y: pd.Series) -> plt.Figure:
    """Cross-validated ROC curves and the AUC null for each candidate signature."""
    from sklearn.metrics import roc_curve

    short = {
        "Whole-compartment composition (11 cell types, CLR)": "composition (11 cell types)",
        "CD8 state composition (CLR)": "CD8 state composition",
        "TCF7+ fraction of CD8 cells [published comparator]": "TCF7+ CD8 fraction (published)",
        "CD8 baseline gene program (nested, 50 genes)": "CD8 baseline gene program",
        "CD8 on-treatment program scored at baseline (nested, 50 genes)":
            "CD8 on-treatment program",
    }
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3),
                             gridspec_kw={"width_ratios": [1.0, 1.0], "wspace": 0.75})
    palette = ["#1f4e79", "#4a90c2", "#D55E00", "#2e7d6f", "#9467bd"]
    ax = axes[0]
    ax.plot([0, 1], [0, 1], ls=(0, (3, 3)), lw=0.9, color="#999999")
    ax.text(0.62, 0.55, "chance", fontsize=6, color="#999999", rotation=37)
    for (name, s), col in zip(scores.items(), palette):
        yy = y.loc[s.index]
        fpr, tpr, _ = roc_curve(yy.to_numpy(int), s.to_numpy())
        auc = float(perf.set_index("signature").loc[name, "auc_loo"])
        ax.plot(fpr, tpr, lw=1.7, color=col,
                label=f"{short.get(name, name)}  (AUC {auc:.2f})")
    ax.set_xlabel("false positive rate"); ax.set_ylabel("true positive rate")
    ax.set_title("Leave-one-patient-out ROC, baseline biopsies")
    ax.legend(frameon=False, fontsize=6.0, loc="upper center",
              bbox_to_anchor=(0.5, -0.17), handlelength=1.4)
    ax.set_aspect("equal"); ax.margins(0.02)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)

    ax2 = axes[1]
    perf2 = perf.sort_values("auc_loo")
    ypos = np.arange(len(perf2))
    cols = {n: c for n, c in zip(scores, palette)}
    for yi, r in zip(ypos, perf2.itertuples()):
        ax2.plot([r.ci95_lo, r.ci95_hi], [yi, yi], lw=1.4,
                 color=cols.get(r.signature, "#444444"), alpha=0.55)
        ax2.scatter([r.auc_loo], [yi], s=26, zorder=3, color=cols.get(r.signature, "#444444"))
        ax2.text(1.02, yi, f"p={r.permutation_p:.3f}", transform=ax2.get_yaxis_transform(),
                 fontsize=6, va="center")
    ax2.axvline(0.5, ls=(0, (3, 3)), lw=0.9, color="#999999")
    ax2.set_yticks(ypos)
    ax2.set_yticklabels([short.get(s, s) for s in perf2["signature"]], fontsize=6.2)
    ax2.set_xlabel("cross-validated AUC (bar = 95% bootstrap CI)")
    ax2.set_xlim(0.0, 1.0)
    ax2.set_title("No signature separates responders beyond the permutation null")
    for s_ in ("top", "right"):
        ax2.spines[s_].set_visible(False)
    return fig


def qc_figure(obs_before: pd.DataFrame, qc: dict) -> plt.Figure:
    """Distributions of the two QC metrics with the applied thresholds."""
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.4))
    specs = [("n_genes", "genes detected per cell", qc["min_genes_per_cell"], "min"),
             ("pct_mito", "% of TPM from mitochondrial genes", qc["max_pct_mito"], "max")]
    for ax, (col, label, thr, kind) in zip(axes, specs):
        v = obs_before[col].to_numpy()
        ax.hist(v, bins=60, color="#4a90c2", lw=0)
        ax.axvline(thr, color="#D55E00", lw=1.3)
        drop = (v < thr).mean() if kind == "min" else (v > thr).mean()
        ax.text(thr, ax.get_ylim()[1] * 0.95,
                f"  {kind} {thr:g}\n  drops {drop*100:.1f}% of cells",
                fontsize=6.5, va="top", color="#D55E00")
        ax.set_xlabel(label); ax.set_ylabel("cells")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.suptitle("QC thresholds remove 11% of cells, almost all on mitochondrial fraction",
                 fontsize=8)
    fig.tight_layout()
    return fig
