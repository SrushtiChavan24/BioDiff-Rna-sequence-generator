"""
Thermodynamic Stability Scoring Service.

Computes a composite Therapeutic Stability Score integrating:
- Thermodynamic stability (MFE per nucleotide)
- GC content (thermal resilience indicator)
- Structural confidence (fold reliability)
- Paired fraction (structural compactness)

The score is normalized to a 0–100 scale where higher values indicate
greater suitability as a therapeutic RNA candidate.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from app.utils.sequence_utils import compute_gc_content


@dataclass
class StabilityResult:
    """Internal stability metrics (before Pydantic conversion)."""

    therapeutic_stability_score: float
    gc_content: float
    mfe_per_nucleotide: float
    ensemble_diversity: float


def _normalize_mfe_score(mfe_per_nt: float) -> float:
    """
    Normalize MFE per nucleotide to a 0–1 score.

    More negative MFE → higher score (more stable).
    Typical range: -0.5 to -0.05 kcal/mol/nt for RNA.
    """
    # Sigmoid-like mapping centered around -0.25 kcal/mol/nt
    score = 1.0 / (1.0 + math.exp(8.0 * (mfe_per_nt + 0.15)))
    return min(1.0, max(0.0, score))


def _gc_content_score(gc: float) -> float:
    """
    Score GC content for therapeutic suitability.

    Optimal GC is 0.40–0.60 for most therapeutic RNAs.
    Penalize extreme values.
    """
    # Gaussian centered at 0.50 with σ = 0.15
    return math.exp(-((gc - 0.50) ** 2) / (2 * 0.15**2))


def _confidence_weight(confidence: float) -> float:
    """
    Weight contribution based on structural confidence.

    High confidence → full weight.
    Low confidence → reduced but non-zero contribution.
    """
    return 0.3 + 0.7 * confidence


def compute_stability_metrics(
    sequence: str,
    mfe: float,
    structural_confidence: float,
    paired_fraction: float,
    ensemble_diversity: float,
) -> StabilityResult:
    """
    Compute comprehensive stability metrics for an RNA sequence.

    Args:
        sequence: Validated RNA sequence.
        mfe: Minimum free energy in kcal/mol.
        structural_confidence: Confidence in the predicted structure (0–1).
        paired_fraction: Fraction of nucleotides in base pairs (0–1).
        ensemble_diversity: Structural ensemble diversity measure.

    Returns:
        StabilityResult with all metrics and composite score.
    """
    n = len(sequence)
    gc = compute_gc_content(sequence)
    mfe_per_nt = mfe / n if n > 0 else 0.0

    # Component scores (each 0–1)
    thermo_score = _normalize_mfe_score(mfe_per_nt)
    gc_score = _gc_content_score(gc)
    compactness_score = min(1.0, paired_fraction / 0.6)  # Normalize to ~60% paired

    # Diversity penalty: lower diversity → more therapeutically predictable
    diversity_penalty = 1.0 - min(1.0, ensemble_diversity)

    # Weighted composite (weights reflect therapeutic relevance)
    confidence_w = _confidence_weight(structural_confidence)
    composite = (
        0.35 * thermo_score
        + 0.25 * gc_score
        + 0.15 * compactness_score
        + 0.15 * diversity_penalty
        + 0.10 * structural_confidence
    ) * confidence_w

    # Scale to 0–100
    therapeutic_stability_score = round(composite * 100, 2)

    return StabilityResult(
        therapeutic_stability_score=min(100.0, max(0.0, therapeutic_stability_score)),
        gc_content=round(gc, 4),
        mfe_per_nucleotide=round(mfe_per_nt, 4),
        ensemble_diversity=round(ensemble_diversity, 4),
    )
