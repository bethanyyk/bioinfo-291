"""
End-to-end demonstration of the statistical-review skill.

Builds a synthetic manuscript PDF and its raw data, with three planted problems
of the kinds this skill exists to find, then runs the full review and writes the
figure spine, the claim ledger and a filled report.

Planted problems:

  C1  pseudoreplication -- 96 measurement fields from 8 mice analysed as 96
      independent observations; the paper reports t(94) while stating n = 4
      per group
  C2  resolution floor -- a paired design with 5 subjects, where the smallest
      two-sided sign-test p is 0.0625 and the paper reports p = 0.040
  C3  no supporting data -- a 12-month claim with only 4-week source data

Run from the skill root:

    python examples/run_demo_review.py --outdir examples/demo_output
"""

import argparse
import hashlib
import os
import sys

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kernel import (audit_reported_statistics, build_review_figures, describe_data,
                    load_paper, review_claim, review_summary_table, verdict_counts)
from utilities.assumption_checks import uncheckable_assumptions_note

RANDOM_STATE = 0
N_PERM = 10000
N_BOOT = 10000

# The synthetic manuscript. Written so the extractors have something realistic to
# parse: reported statistics with df, an unlabelled error bar, a design whose
# stated n contradicts its reported degrees of freedom, and enough p-values for
# the multiplicity diagnostics to run.
MANUSCRIPT_PARAGRAPHS = [
    "Compound X reduces tumour burden in a syngeneic model",
    "",
    "Abstract. Treatment with compound X reduced tumour volume relative to vehicle "
    "control (t(94) = 3.41, p = 0.001). The effect was reproduced in a paired "
    "crossover cohort (p = 0.040) and persisted at 12 months.",
    "",
    "Methods. Female mice were housed in the institutional vivarium. Animals were "
    "assigned to treatment or vehicle control, n = 4 per group. Tumour volume was "
    "measured in 12 independent fields per animal. Data are shown as mean +/- error "
    "from n = 4 independent experiments. Groups were compared with Student's t-test. "
    "Paired measurements in the crossover cohort were compared with a Wilcoxon "
    "signed-rank test. No correction for multiple comparisons was applied.",
    "",
    "Results. Tumour volume was lower in treated animals than in vehicle controls "
    "(t(94) = 3.41, p = 0.001). Marker abundance correlated with final volume "
    "(r(94) = 0.26, p = 0.011). Across the marker panel, differences were observed "
    "for CD8 (p = 0.012), CD4 (p = 0.047), FOXP3 (p = 0.049), GZMB (p = 0.046), "
    "PRF1 (p = 0.048), IFNG (p = 0.043), TNF (p = 0.049), IL2 (p = 0.044), "
    "IL10 (p = 0.21), TGFB1 (p = 0.33), VEGFA (p = 0.52), PDCD1 (p = 0.61), "
    "CTLA4 (p = 0.74), LAG3 (p = 0.81), HAVCR2 (p = 0.09), TIGIT (p = 0.15), "
    "KLRG1 (p = 0.28), TBX21 (p = 0.37), EOMES (p = 0.55) and PRDM1 (p = 0.68). "
    "In the paired crossover cohort, volume fell in 4 of 5 subjects (p = 0.040). "
    "Counts of responding animals differed between arms (chi2(1) = 2.10, p = 0.012).",
    "",
    "Discussion. These data indicate that compound X is efficacious across the immune "
    "compartment and that the effect is durable to 12 months.",
]


def _wrap(text: str, width: int):
    words, lines, cur = text.split(), [], ""
    for w in words:
        if len(cur) + len(w) + 1 > width:
            lines.append(cur); cur = w
        else:
            cur = f"{cur} {w}".strip()
    if cur:
        lines.append(cur)
    return lines


def build_demo_pdf(path: str) -> None:
    """Render the synthetic manuscript to a PDF with a real text layer."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with PdfPages(path) as pdf:
        fig = plt.figure(figsize=(8.27, 11.69))
        y = 0.96
        for para in MANUSCRIPT_PARAGRAPHS:
            if not para:
                y -= 0.015
                continue
            for line in _wrap(para, 92):
                fig.text(0.08, y, line, fontsize=9, va="top")
                y -= 0.018
            y -= 0.012
        pdf.savefig(fig)
        plt.close(fig)


def build_demo_data(rng: np.random.Generator):
    """
    Two source datasets matching two of the paper's figures.

    `fields` is the pseudoreplicated one: 12 measurement fields per animal from
    8 animals, with between-animal variance large relative to within-animal, so
    a row-level test is far more confident than the design allows.

    `crossover` is the 5-subject paired cohort.
    """
    rows = []
    for unit in range(8):
        arm = "compound X" if unit < 4 else "vehicle"
        animal_effect = rng.normal(-0.35 if arm == "compound X" else 0.0, 0.55)
        for field in range(12):
            rows.append({
                "mouse": f"M{unit + 1:02d}",
                "arm": arm,
                "field": field + 1,
                "volume_mm3": float(np.exp(5.2 + animal_effect + rng.normal(0, 0.12))),
            })
    fields = pd.DataFrame(rows)

    subjects = [f"S{i + 1}" for i in range(5)]
    baseline = rng.normal(180, 25, 5)
    change = np.array([-18.0, -12.0, -9.0, -15.0, 4.0])  # 4 of 5 improve
    crossover = pd.DataFrame({
        "subject": np.repeat(subjects, 2),
        "arm": ["baseline", "on compound X"] * 5,
        "volume_mm3": np.ravel([[b, b + c] for b, c in zip(baseline, change)]),
    })
    return fields, crossover


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def main(outdir: str) -> dict:
    rng = np.random.default_rng(RANDOM_STATE)
    raw_dir = os.path.join(outdir, "data")
    fig_dir = os.path.join(outdir, "figures")
    res_dir = os.path.join(outdir, "results")
    for d in (raw_dir, fig_dir, res_dir):
        os.makedirs(d, exist_ok=True)

    pdf_path = os.path.join(raw_dir, "manuscript.pdf")
    build_demo_pdf(pdf_path)
    fields, crossover = build_demo_data(rng)
    fields_path = os.path.join(raw_dir, "figure1_fields.csv")
    cross_path = os.path.join(raw_dir, "figure2_crossover.csv")
    fields.to_csv(fields_path, index=False)
    crossover.to_csv(cross_path, index=False)

    # ---- stage 0: establish what we were given -----------------------------
    paper = load_paper(pdf_path)

    # ---- stages 2 and 6, from the PDF alone --------------------------------
    audit = audit_reported_statistics(paper["inventory"])

    # ---- stages 1-5: review each claim -------------------------------------
    c1 = review_claim(
        fields, value_col="volume_mm3", group_col="arm", unit_col="mouse",
        claim_id="C1", claim_short="compound X reduces tumour volume",
        p_reported=0.001, n_perm=N_PERM, n_boot=N_BOOT, random_state=RANDOM_STATE,
        figure_path=os.path.join(fig_dir, "panel_C1.png"),
    )
    # the same claim at the unit level, which is the n the design supports
    per_mouse = fields.groupby(["mouse", "arm"], as_index=False)["volume_mm3"].mean()
    c1_unit = review_claim(
        per_mouse, value_col="volume_mm3", group_col="arm",
        claim_id="C1u", claim_short="C1 re-tested at the animal level",
        p_reported=0.001, n_perm=N_PERM, n_boot=N_BOOT, random_state=RANDOM_STATE,
        figure_path=os.path.join(fig_dir, "panel_C1u.png"),
    )
    c2 = review_claim(
        crossover, value_col="volume_mm3", group_col="arm", paired=True,
        pair_col="subject", claim_id="C2",
        claim_short="effect reproduced in paired crossover",
        p_reported=0.040, n_perm=N_PERM, n_boot=N_BOOT, random_state=RANDOM_STATE,
        figure_path=os.path.join(fig_dir, "panel_C2.png"),
    )
    claims = [c1, c1_unit, c2]

    # ---- the fixed figure spine --------------------------------------------
    figs = build_review_figures(
        fields, "mouse", "arm", "volume_mm3", claims=claims, audit=audit,
        outdir=fig_dir, value_label="tumour volume (mm3)",
        reported_error_label="unstated ('mean +/- error')",
    )

    # ---- ledger, including the claim with no data --------------------------
    ledger = review_summary_table(claims)
    ledger = pd.concat([ledger, pd.DataFrame([{
        "claim_id": "C3",
        "claim_short": "effect persists at 12 months",
        "verdict": "not-assessable",
        "effect_reported": np.nan, "effect_recomputed": np.nan,
        "ci_low": np.nan, "ci_high": np.nan,
        "p_reported": np.nan, "p_recomputed": np.nan,
        "reasons": "no longitudinal data supplied; source data cover the 4-week endpoint only",
    }])], ignore_index=True)
    ledger.to_csv(os.path.join(res_dir, "claim_ledger.csv"), index=False)
    audit["statcheck"].to_csv(os.path.join(res_dir, "statcheck.csv"), index=False)
    if audit["multiplicity"]:
        audit["multiplicity"]["summary"].to_csv(
            os.path.join(res_dir, "multiplicity_summary.csv"), index=False)

    # rebuild the ledger figure so the not-assessable claim appears in it
    from utilities.review_figures import fig_claim_ledger
    fig_claim_ledger(ledger, path=figs["F4_claim_ledger"],
                     effect_label="difference in tumour volume (mm3)")

    desc = describe_data(fields, "volume_mm3", "arm", unit_col="mouse")
    report = write_report(outdir, paper, audit, desc, claims, ledger, figs,
                          {"fields": fields_path, "crossover": cross_path, "pdf": pdf_path})

    return {"paper": paper, "audit": audit, "claims": claims, "ledger": ledger,
            "figures": figs, "report": report, "describe": desc}


def write_report(outdir, paper, audit, desc, claims, ledger, figs, paths) -> str:
    """Fill the report template with this review's numbers."""
    ind = desc["independence"]
    mult = audit["multiplicity"]
    jb = audit["just_below"]
    c1, c1u, c2 = claims
    counts = ledger["verdict"].value_counts()

    def rel(p):
        return os.path.relpath(p, outdir).replace("\\", "/")

    lines = [
        "# Statistical Review - Compound X reduces tumour burden in a syngeneic model",
        "",
        "*Reviewer: statistical-review skill (demonstration run) · "
        "Protocol: statistical-review v1*",
        "", "---", "",
        "## 1. Inputs and coverage",
        "",
        "| Input | Value |",
        "|-------|-------|",
        f"| Manuscript | `{os.path.basename(paths['pdf'])}` - "
        f"{paper['text_quality']['n_pages']} page(s), text layer "
        f"{'usable' if paper['text_quality']['usable_text_layer'] else 'NOT usable'} "
        f"({paper['text_quality']['chars_per_page']:.0f} chars/page) |",
        f"| Raw data (Figure 1) | `{os.path.basename(paths['fields'])}` - "
        f"{desc['n_rows']} rows, sha256 `{sha256(paths['fields'])[:16]}...` |",
        f"| Raw data (Figure 2) | `{os.path.basename(paths['crossover'])}` - "
        f"sha256 `{sha256(paths['crossover'])[:16]}...` |",
        f"| Parsed from the PDF | {len(paper['inventory']['tests'])} test statistics, "
        f"{len(paper['inventory']['p_values'])} p-values |",
        "| Claims with no supporting data | C3 - reported as `not-assessable` |",
        f"| Seeds | `random_state={RANDOM_STATE}`, `n_perm={N_PERM}`, `n_boot={N_BOOT}` |",
        "",
        "This review covers 3 of the paper's claims; C3 could not be assessed with the "
        "supplied data.",
        "", "---", "",
        "## 2. Summary of findings",
        "",
        "| Verdict | Claims |",
        "|---------|--------|",
    ]
    for verdict, n in counts.items():
        lines.append(f"| `{verdict}` | {n} |")
    lines += [
        "",
        f"1. **C1 is pseudoreplicated.** The paper states *n = 4 per group* but reports "
        f"*t(94)*, degrees of freedom consistent with a test on all "
        f"{ind['n_observations']} measurement fields. Fields cluster strongly within animal "
        f"(ICC = {ind['icc']:.2f}, {ind['obs_per_unit_mean']:.0f} per animal), giving an "
        f"effective n of about {ind['effective_n']:.1f}. Re-tested at the animal level (C1u) "
        f"the permutation p is {c1u['p_recomputed']:.3f} against the reported 0.001.",
        f"2. **C2 reports a p-value its design cannot produce.** With "
        f"{c2['permutation']['n1']} pairs the sign-flip permutation floor is p >= "
        f"{c2['permutation']['min_achievable_p']:.4f} and the sign test cannot go below "
        f"0.0625; the paper reports p = 0.040.",
        f"3. **No correction was applied across {mult['m_tests']} reported p-values**, where "
        f"about {mult['expected_false_positives_uncorrected']:.1f} significant results are "
        f"expected by chance alone. {int(mult['summary'].loc[0, 'n_significant'])} uncorrected "
        f"hits become {int(mult['summary'].loc[2, 'n_significant'])} under Benjamini-Hochberg "
        f"and {int(mult['summary'].loc[1, 'n_significant'])} under Bonferroni.",
        "", "---", "",
        "## 3. Claim ledger",
        "",
        f"![F4 - claim ledger]({rel(figs['F4_claim_ledger'])})",
        "",
        "*Recomputed effect (filled) against reported (open) with 95% bootstrap CI; right, "
        "reported against recomputed p-values on log axes with alpha = 0.05 marked.*",
        "",
        "| id | claim | verdict | effect recomp. | 95% CI | p rep. | p recomp. |",
        "|----|-------|---------|----------------|--------|--------|-----------|",
    ]
    for r in ledger.itertuples():
        ci = f"[{r.ci_low:.1f}, {r.ci_high:.1f}]" if pd.notna(r.ci_low) else "-"
        lines.append(
            f"| {r.claim_id} | {r.claim_short} | `{r.verdict}` | "
            f"{'-' if pd.isna(r.effect_recomputed) else format(r.effect_recomputed, '.1f')} | "
            f"{ci} | "
            f"{'-' if pd.isna(r.p_reported) else format(r.p_reported, '.3f')} | "
            f"{'-' if pd.isna(r.p_recomputed) else format(r.p_recomputed, '.4f')} |"
        )
    lines += ["", "Reasons, per claim:", ""]
    for r in ledger.itertuples():
        lines.append(f"- **{r.claim_id}**: {r.reasons}")
    lines += [
        "", "---", "",
        "## 4. Design and sampling - stages 1-2",
        "",
        f"![F1 - design map]({rel(figs['F1_design_map'])})",
        "",
        "*Left: measurement fields per animal. Right: independent animals per arm, the n "
        "available to a test.*",
        "",
        f"![F2 - raw data]({rel(figs['F2_raw_behind_claim'])})",
        "",
        "*Every field with its box summary; right panel shows mean +/- SD against "
        "mean +/- SEM on the same axis.*",
        "",
        "**Population.** The methods state only that animals were female and vivarium-housed. "
        "Strain, age and facility are not given, so the population the claim covers cannot be "
        "pinned down from the manuscript.",
        "",
        f"**Unit of replication.** {ind['note']} The reported *t(94)* is consistent with a "
        f"test on measurement fields; a two-sample test on the {ind['n_units']} animals has "
        f"{ind['n_units'] - 2} degrees of freedom.",
        "",
        "**Uncertainty reporting.** The paper writes \"mean +/- error\" without saying which "
        "quantity the bars represent. SD and SEM differ by sqrt(n) here, so the figure cannot "
        "be read as drawn.",
        "", "---", "",
        "## 5. Method and assumptions - stages 3-4",
        "",
        f"![F3 - distribution shape]({rel(figs['F3_assumption_shape'])})",
        "",
        "*Q-Q per arm with n, skewness, excess kurtosis and the Shapiro-Wilk p-value. "
        "Descriptive, not a pass/fail gate.*",
        "",
        f"**Test selection.** The paper names Student's equal-variance t-test. "
        f"{desc['variance']['note']} The crossover cohort is analysed with Wilcoxon "
        "signed-rank, which assumes the paired differences are symmetric; the sign test makes "
        "no such assumption and is reported alongside below.",
        "",
        "**Three routes, C1 at the animal level (C1u):**",
        "",
        "| route | test | p | assumes |",
        "|-------|------|---|---------|",
    ]
    for r in c1u["routes"].itertuples():
        lines.append(f"| {r.route} | {r.test} | {r.p_value:.4f} | {r.assumes} |")
    lines += [
        "",
        "**Three routes, C2 (paired, 5 subjects):**",
        "",
        "| route | test | p | smallest p this design allows |",
        "|-------|------|---|-------------------------------|",
    ]
    for r in c2["routes"].itertuples():
        lines.append(f"| {r.route} | {r.test} | {r.p_value:.4f} | {r.min_achievable_p:.4f} |")
    lines += [
        "", "---", "",
        "## 6. Resolution and evidence - stage 5",
        "",
        f"![F5 - null distribution]({rel(figs['F5_null_distribution'])})",
        "",
        "*Permutation null with the observed statistic marked, annotated with the p-value, "
        "the resolution limit and the design floor where it binds; right, the bootstrap "
        "distribution of the estimate with its interval.*",
        "",
        f"**C1 (field level).** {c1['permutation']['note']}",
        "",
        f"**C1u (animal level).** 4 vs 4 animals gives a two-sided permutation floor of "
        f"{c1u['permutation']['min_achievable_p']:.4f}; observed p = "
        f"{c1u['p_recomputed']:.4f}. Bootstrap difference in means "
        f"{c1u['bootstrap']['point_estimate']:.1f} mm3, 95% CI "
        f"[{c1u['bootstrap']['ci_low']:.1f}, {c1u['bootstrap']['ci_high']:.1f}].",
        "",
        f"**C2 (paired, 5 subjects).** {c2['permutation']['note']}",
        "",
        "**Shuffling scheme.** The manuscript does not state how, or whether, any permutation "
        "was performed; the p-values above come from this review, permuting labels across "
        "animals for C1u and flipping signs within subject for C2.",
        "", "---", "",
        "## 7. Multiplicity - stage 6",
        "",
        f"![F6 - multiplicity]({rel(figs['F6_multiplicity'])})",
        "",
        "*p-value histogram against the uniform expectation with the just-below-threshold "
        "window highlighted; right, hits surviving each correction against the number "
        "expected by chance.*",
        "",
        "| correction | controls | n significant | threshold on raw p |",
        "|------------|----------|---------------|--------------------|",
    ]
    for r in mult["summary"].itertuples():
        lines.append(f"| {r.correction} | {r.controls} | {r.n_significant} | "
                     f"{r.threshold_on_raw_p:.4g} |")
    lines += [
        "",
        f"**Counting.** {mult['note']} The methods state that no correction was applied.",
        "",
        f"**Distribution of reported p-values.** {jb['note']} This diagnostic describes the "
        "set of reported p-values and says nothing about any individual result; it cannot "
        "distinguish analytic flexibility from legitimate exploratory analysis, and rounding "
        "conventions alone can produce apparent pile-up. The full set of comparisons run, or a "
        "pre-specified analysis plan, would resolve it.",
        "", "---", "",
        "## 8. Per-claim supplement",
        "",
    ]
    for claim in claims:
        lines += [
            f"![{claim['claim_id']}]({rel(claim['figure_path'])})",
            "",
            f"*{claim['claim_id']} - {claim['claim_short']}; verdict `{claim['verdict']}`. "
            "Raw observations, permutation null, and the three routes. Shaded region marks "
            "p-values unreachable for this design.*",
            "",
        ]
    lines += [
        "---", "",
        "## 9. What was not assessed",
        "",
        uncheckable_assumptions_note(),
        "",
        "Additionally:",
        "",
        "- C3 (\"effect persists at 12 months\"): no longitudinal data supplied.",
        "- The paper's own model fit could not be reconstructed; only the group comparisons "
        "were re-derived.",
        "- Blinding and randomization are not described in the methods, so neither was "
        "assessed.",
        "", "---", "",
        "## 10. Reproduction",
        "",
        "```bash",
        "python examples/run_demo_review.py --outdir examples/demo_output",
        "```",
        "",
        "| | |",
        "|--|--|",
        f"| Seeds | `random_state={RANDOM_STATE}` |",
        f"| Permutations / resamples | `n_perm={N_PERM}`, `n_boot={N_BOOT}` |",
        f"| Package versions | numpy {np.__version__}, pandas {pd.__version__} |",
        f"| Manuscript checksum | sha256 `{sha256(paths['pdf'])[:16]}...` |",
        "| Skill version | statistical-review v1 |",
        "", "---", "",
        "*Findings in this report are arithmetic statements about reported and recomputed "
        "values. Where they differ, the cause is often an input the review could not see; the "
        "inputs used are named above so that discrepancies can be resolved.*",
        "",
    ]

    path = os.path.join(outdir, "REVIEW_REPORT.md")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    return path


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(os.path.dirname(__file__), "demo_output"))
    args = ap.parse_args()
    out = main(args.outdir)
    print(f"report:  {out['report']}")
    print(f"figures: {len(out['figures'])}")
    print(out["ledger"][["claim_id", "verdict", "p_reported", "p_recomputed"]]
          .to_string(index=False))
