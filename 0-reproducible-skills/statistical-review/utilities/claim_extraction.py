"""
Extract statistical claims from a paper PDF.

A "claim" is any sentence in which the paper asserts something about a
population on the basis of a sample: a reported test, a p-value, an interval,
an effect size, or a sample size.

This module only *locates and parses* claims. It makes no judgement; judgement
happens in recompute.py (does the arithmetic hold?) and in the review protocol
(does the claim follow from the data?).

Nothing here is a substitute for reading the paper. Regexes miss claims stated
in prose ("treatment clearly improved survival") and claims that live only in
figures. Treat the output as a starting inventory to be completed by hand.
"""

import re
from typing import List, Optional, Sequence

import numpy as np
import pandas as pd

# Statistics that can be re-derived into a p-value, mapped to their df arity.
TEST_STAT_ARITY = {
    "t": 1,
    "F": 2,
    "r": 1,
    "chi2": 1,
    "z": 0,
}

# Keywords that signal how uncertainty was reported. Which one a paper uses
# decides whether an error bar is readable at all.
ERROR_BAR_KEYWORDS = (
    "standard error",
    "s.e.m",
    "sem",
    "standard deviation",
    "s.d.",
    "error bar",
    "confidence interval",
    "credible interval",
    "interquartile",
    "iqr",
    "whisker",
    "±",
)

# Keywords that signal the unit of analysis, i.e. what n actually counts.
REPLICATION_KEYWORDS = (
    "biological replicate",
    "technical replicate",
    "independent experiment",
    "independent repeat",
    "per group",
    "per animal",
    "per mouse",
    "per patient",
    "per donor",
    "per well",
    "pooled",
    "cells from",
)

# Keywords that signal a multiplicity-relevant analysis.
MULTIPLICITY_KEYWORDS = (
    "bonferroni",
    "holm",
    "benjamini",
    "false discovery",
    "fdr",
    "q-value",
    "q value",
    "adjusted p",
    "corrected for multiple",
    "multiple comparison",
    "multiple testing",
    "post hoc",
    "post-hoc",
)

_SENTENCE_SPLIT = re.compile(r"(?<=[.!?])\s+(?=[A-Z(\[])")

# t(18) = 2.31 / F(2, 27) = 4.10 / r(48) = .31 / chi2(1) = 6.2 / Z = 1.96
_STAT_RE = re.compile(
    r"(?P<stat>\bt\b|\bF\b|\br\b|\bZ\b|\bz\b|chi2|chi-square[d]?|X2)"
    r"\s*(?:\(\s*(?P<df1>\d+(?:[.,]\d+)?)\s*(?:[,;]\s*(?P<df2>\d+(?:[.,]\d+)?)\s*)?\))?"
    r"\s*(?P<rel>[=<>~]{1,2})\s*"
    r"(?P<value>[-+]?\d*\.?\d+)",
)

# df reported separately: "t = 2.31, df = 18"
_DF_RE = re.compile(
    r"\bd\.?f\.?\s*(?:[=:]|\s)\s*(?P<df1>\d+(?:\.\d+)?)(?:\s*[,;]\s*(?P<df2>\d+(?:\.\d+)?))?"
)

# p = .03 / p < 0.001 / P-value = 1.2e-5 / p <= 0.05
_P_RE = re.compile(
    r"\bp[\s\-]?(?:value|val)?\s*(?P<rel>[=<>]{1,2}|is|was)\s*"
    r"(?P<value>\d*\.?\d+(?:\s*[xX]\s*10\s*-?\s*\d+|[eE][-+]?\d+)?)",
    re.IGNORECASE,
)

# n = 12 / N = 3 mice / (n = 24 patients)
_N_RE = re.compile(r"\bn\s*=\s*(?P<n>\d+)\s*(?P<unit>[a-zA-Z][a-zA-Z\- ]{0,24})?", re.IGNORECASE)

# 95% CI [0.2, 1.8] / 95% confidence interval of 0.2 to 1.8
_CI_RE = re.compile(
    r"(?P<level>\d{2}(?:\.\d+)?)\s*%\s*(?:CI|confidence interval|credible interval)"
    r"\s*(?:of|was|:|=)?\s*[\[\(]?\s*"
    r"(?P<lo>[-+]?\d*\.?\d+)\s*(?:,|to|-)\s*(?P<hi>[-+]?\d*\.?\d+)\s*[\]\)]?",
    re.IGNORECASE,
)

# mean +/- error, where the error is only interpretable if the paper says what it is
_PM_RE = re.compile(r"(?P<mean>[-+]?\d*\.?\d+)\s*(?:±|\+/-|\+-)\s*(?P<err>\d*\.?\d+)")

# Named tests, so the review can check that the test matches the design.
TEST_NAME_PATTERNS = {
    "student_t": r"student'?s?\s+t[-\s]?test|unpaired t[-\s]?test|two[-\s]sample t[-\s]?test",
    "welch_t": r"welch'?s?|unequal variance|heteroscedastic t",
    "paired_t": r"paired t[-\s]?test",
    "one_sample_t": r"one[-\s]sample t[-\s]?test",
    "anova": r"\banova\b|analysis of variance",
    "mann_whitney": r"mann[-\s]?whitney|wilcoxon rank[-\s]?sum",
    "wilcoxon_signed_rank": r"wilcoxon signed[-\s]?rank|signed[-\s]?rank test",
    "sign_test": r"\bsign test\b",
    "kruskal": r"kruskal[-\s]?wallis",
    "permutation": r"permutation test|randomi[sz]ation test|label shuffl",
    "bootstrap": r"bootstrap",
    "chi2": r"chi[-\s]?square|chi2 test|fisher'?s exact",
    "correlation": r"pearson|spearman|kendall",
    "regression": r"linear model|regression|mixed[-\s]?effect|generali[sz]ed linear",
    "survival": r"log[-\s]?rank|cox proportional|kaplan[-\s]?meier",
}

_UNICODE_FIX = {
    "\u2212": "-", "\u2013": "-", "\u2014": "-", "\u2010": "-", "\u00a0": " ",
    "\u03c7\u00b2": "chi2", "\u03c7\u0032": "chi2", "\u2264": "<=", "\u2265": ">=",
    "\u00d7": "x",
}


def normalize_text(text: str) -> str:
    """Normalize unicode minus signs, non-breaking spaces and chi-square glyphs."""
    for bad, good in _UNICODE_FIX.items():
        text = text.replace(bad, good)
    return re.sub(r"[ \t]+", " ", text)


def read_pdf_text(pdf_path: str) -> List[str]:
    """
    Extract the text layer of a PDF, one string per page.

    Uses pypdfium2. Pages with no text layer come back as empty strings -- a
    scanned PDF yields all-empty pages, which the caller must report as
    "claims not extractable" rather than "no claims found".
    """
    import pypdfium2 as pdfium

    pdf = pdfium.PdfDocument(pdf_path)
    pages = []
    try:
        for i in range(len(pdf)):
            page = pdf[i]
            textpage = page.get_textpage()
            pages.append(normalize_text(textpage.get_text_range() or ""))
            textpage.close()
            page.close()
    finally:
        pdf.close()
    return pages


def pdf_text_quality(pages: Sequence[str]) -> dict:
    """
    Report whether the PDF has a usable text layer.

    Call this before trusting an empty claim table. A scanned or image-only PDF
    is a finding ("claims could not be extracted"), not an absence of claims.
    """
    n_chars = sum(len(p) for p in pages)
    empty = sum(1 for p in pages if len(p.strip()) < 20)
    return {
        "n_pages": len(pages),
        "n_chars": n_chars,
        "empty_pages": empty,
        "chars_per_page": n_chars / max(len(pages), 1),
        "usable_text_layer": bool(pages) and (n_chars / max(len(pages), 1)) > 200,
    }


def split_sentences(text: str) -> List[str]:
    """Split page text into sentence-ish units for context capture."""
    text = re.sub(r"\s*\n\s*", " ", text)
    return [s.strip() for s in _SENTENCE_SPLIT.split(text) if s.strip()]


def _to_float(value: Optional[str]) -> float:
    if value is None:
        return np.nan
    v = value.replace("\u2212", "-").replace(",", ".").strip()
    v = re.sub(r"\s*[xX]\s*10\s*-\s*(\d+)", r"e-\1", v)
    v = re.sub(r"\s*[xX]\s*10\s*(\d+)", r"e\1", v)
    try:
        return float(v)
    except ValueError:
        return np.nan


def _canon_stat(raw: str) -> str:
    s = raw.strip().lower()
    if s.startswith("chi") or s == "x2":
        return "chi2"
    return {"t": "t", "f": "F", "r": "r", "z": "z"}.get(s, s)


def find_reported_tests(text: str, page: Optional[int] = None) -> pd.DataFrame:
    """
    Find reported test statistics and the p-value attached to each.

    Returns one row per statistic found, with columns: page, sentence,
    stat_type, stat_value, df1, df2, p_relation, p_reported, tail_hint.

    This is the input to recompute.statcheck_table(), which checks whether the
    reported p-value follows from the reported statistic and df.
    """
    rows = []
    for sent in split_sentences(text):
        stats_found = list(_STAT_RE.finditer(sent))
        if not stats_found:
            continue
        ps = list(_P_RE.finditer(sent))
        df_extra = _DF_RE.search(sent)
        for sm in stats_found:
            stat_type = _canon_stat(sm.group("stat"))
            if stat_type not in TEST_STAT_ARITY:
                continue
            df1 = _to_float(sm.group("df1"))
            df2 = _to_float(sm.group("df2"))
            if np.isnan(df1) and df_extra is not None:
                df1 = _to_float(df_extra.group("df1"))
                df2 = _to_float(df_extra.group("df2"))
            # attach the nearest p-value that appears after the statistic
            p_match = next((p for p in ps if p.start() > sm.end()), ps[0] if ps else None)
            rows.append({
                "page": page,
                "sentence": sent[:400],
                "stat_type": stat_type,
                "stat_value": _to_float(sm.group("value")),
                "df1": df1,
                "df2": df2,
                "p_relation": p_match.group("rel").lower() if p_match else None,
                "p_reported": _to_float(p_match.group("value")) if p_match else np.nan,
                "tail_hint": 1 if re.search(r"one[-\s]?(tailed|sided)", sent, re.I) else 2,
            })
    return pd.DataFrame(rows, columns=[
        "page", "sentence", "stat_type", "stat_value", "df1", "df2",
        "p_relation", "p_reported", "tail_hint",
    ])


def find_p_values(text: str, page: Optional[int] = None) -> pd.DataFrame:
    """Find every reported p-value, with or without an accompanying statistic."""
    rows = []
    for sent in split_sentences(text):
        for m in _P_RE.finditer(sent):
            rows.append({
                "page": page,
                "sentence": sent[:400],
                "p_relation": m.group("rel").lower(),
                "p_reported": _to_float(m.group("value")),
                "is_bound": m.group("rel") in ("<", ">", "<=", ">="),
            })
    return pd.DataFrame(rows, columns=["page", "sentence", "p_relation", "p_reported", "is_bound"])


def find_sample_sizes(text: str, page: Optional[int] = None) -> pd.DataFrame:
    """
    Find reported sample sizes and whatever noun follows them.

    The noun matters more than the number: "n = 3 independent experiments" and
    "n = 3000 cells from 3 mice" can describe the same experiment, and only the
    first states the unit the test should treat as independent.
    """
    rows = []
    for sent in split_sentences(text):
        for m in _N_RE.finditer(sent):
            unit = (m.group("unit") or "").strip().lower()
            rows.append({
                "page": page,
                "sentence": sent[:400],
                "n": int(m.group("n")),
                "unit_stated": unit or None,
                "replication_language": [k for k in REPLICATION_KEYWORDS if k in sent.lower()] or None,
            })
    return pd.DataFrame(rows, columns=["page", "sentence", "n", "unit_stated", "replication_language"])


def find_intervals(text: str, page: Optional[int] = None) -> pd.DataFrame:
    """Find reported confidence/credible intervals and bare mean +/- error pairs."""
    rows = []
    for sent in split_sentences(text):
        for m in _CI_RE.finditer(sent):
            rows.append({
                "page": page, "sentence": sent[:400], "kind": "interval",
                "level": _to_float(m.group("level")),
                "lo": _to_float(m.group("lo")), "hi": _to_float(m.group("hi")),
                "center": np.nan, "error": np.nan,
                "error_named": [k for k in ERROR_BAR_KEYWORDS if k in sent.lower()] or None,
            })
        for m in _PM_RE.finditer(sent):
            rows.append({
                "page": page, "sentence": sent[:400], "kind": "mean_pm_error",
                "level": np.nan, "lo": np.nan, "hi": np.nan,
                "center": _to_float(m.group("mean")), "error": _to_float(m.group("err")),
                "error_named": [k for k in ERROR_BAR_KEYWORDS if k in sent.lower()] or None,
            })
    return pd.DataFrame(rows, columns=[
        "page", "sentence", "kind", "level", "lo", "hi", "center", "error", "error_named",
    ])


def find_named_tests(text: str, page: Optional[int] = None) -> pd.DataFrame:
    """Find which statistical tests the paper names, and where."""
    rows = []
    low = text.lower()
    for name, pattern in TEST_NAME_PATTERNS.items():
        for m in re.finditer(pattern, low, re.IGNORECASE):
            start = max(0, m.start() - 160)
            rows.append({
                "page": page,
                "test_named": name,
                "matched_text": text[m.start():m.end()],
                "context": normalize_text(text[start:m.end() + 160]).strip()[:400],
            })
    return pd.DataFrame(rows, columns=["page", "test_named", "matched_text", "context"])


def find_keyword_sentences(text: str, keywords: Sequence[str],
                           page: Optional[int] = None) -> pd.DataFrame:
    """Find sentences containing any of `keywords` (error bars, multiplicity, design)."""
    rows = []
    for sent in split_sentences(text):
        hits = [k for k in keywords if k in sent.lower()]
        if hits:
            rows.append({"page": page, "sentence": sent[:400], "keywords": hits})
    return pd.DataFrame(rows, columns=["page", "sentence", "keywords"])


def build_claim_inventory(pages: Sequence[str]) -> dict:
    """
    Run every extractor over every page and return a dict of DataFrames.

    Keys: tests, p_values, sample_sizes, intervals, named_tests, error_bars,
    multiplicity, replication, text_quality.

    Each table carries `page` so the reviewer can go back to the source. The
    inventory is deliberately over-inclusive: it is cheaper to discard a false
    hit than to miss a claim.
    """
    def _concat(fn):
        frames = [fn(p, i + 1) for i, p in enumerate(pages)]
        frames = [f for f in frames if len(f)]
        return pd.concat(frames, ignore_index=True) if frames else fn("", None)

    inv = {
        "tests": _concat(find_reported_tests),
        "p_values": _concat(find_p_values),
        "sample_sizes": _concat(find_sample_sizes),
        "intervals": _concat(find_intervals),
        "named_tests": _concat(find_named_tests),
        "error_bars": _concat(lambda t, page=None: find_keyword_sentences(t, ERROR_BAR_KEYWORDS, page)),
        "multiplicity": _concat(lambda t, page=None: find_keyword_sentences(t, MULTIPLICITY_KEYWORDS, page)),
        "replication": _concat(lambda t, page=None: find_keyword_sentences(t, REPLICATION_KEYWORDS, page)),
    }
    inv["text_quality"] = pd.DataFrame([pdf_text_quality(pages)])
    for key, frame in inv.items():
        if key != "text_quality" and len(frame):
            frame.insert(0, "claim_id", [f"{key[:3].upper()}{i + 1:03d}" for i in range(len(frame))])
    return inv


def inventory_summary(inventory: dict) -> pd.DataFrame:
    """One row per inventory table with its row count, for the report header."""
    return pd.DataFrame(
        [{"table": k, "n_rows": len(v)} for k, v in inventory.items() if k != "text_quality"]
    ).sort_values("n_rows", ascending=False, ignore_index=True)
