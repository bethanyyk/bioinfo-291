"""
Statistical-review utilities.

    claim_extraction  -- locate and parse statistical claims in the paper PDF
    recompute         -- re-derive reported p-values, means and effects
    assumption_checks -- describe distribution shape, spread and clustering
    resampling        -- permutation, sign and bootstrap, with resolution floors
    multiplicity      -- correction recounts and p-value distribution diagnostics
    review_figures    -- the six standard review figures and the per-claim panel

Import the submodules you need:

    from utilities import claim_extraction, recompute, resampling
"""

from . import (
    assumption_checks,
    claim_extraction,
    multiplicity,
    recompute,
    resampling,
    review_figures,
)

__all__ = [
    "claim_extraction",
    "recompute",
    "assumption_checks",
    "resampling",
    "multiplicity",
    "review_figures",
]
