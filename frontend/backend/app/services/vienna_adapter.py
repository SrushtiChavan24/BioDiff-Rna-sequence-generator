from __future__ import annotations

"""
ViennaRNA Integration Adapter.

Attempts to import ViennaRNA (RNAlib) for secondary structure prediction.
If unavailable, falls back to the pure-Python Nussinov engine.

Both paths return a unified StructurePredictionResult.
"""

import logging
from dataclasses import dataclass

logger = logging.getLogger(__name__)

# Attempt ViennaRNA import
_VIENNA_AVAILABLE = False
try:
    import RNA  # ViennaRNA Python bindings
    _VIENNA_AVAILABLE = True
    logger.info("ViennaRNA detected — using RNA.fold() for structure prediction")
except ImportError:
    logger.info("ViennaRNA not installed — using pure-Python Nussinov engine")


@dataclass
class FoldResult:
    dot_bracket: str
    minimum_free_energy: float
    base_pairs: list[list[int]]
    structural_confidence: float
    paired_fraction: float
    ensemble_diversity: float


def fold_with_vienna(sequence: str) -> FoldResult | None:
    """Predict secondary structure using ViennaRNA.

    Returns None if ViennaRNA is unavailable or prediction fails.
    """
    if not _VIENNA_AVAILABLE:
        return None

    try:
        n = len(sequence)
        fc = RNA.fold_compound(sequence)
        mfe_str, mfe = fc.mfe()

        # Parse MFE string into dot-bracket
        dot_bracket = mfe_str.strip()

        # Extract base pairs from dot-bracket using a stack
        stack: list[int] = []
        pairs: list[list[int]] = []
        for i, ch in enumerate(dot_bracket):
            if ch == "(":
                stack.append(i)
            elif ch == ")":
                if stack:
                    j = stack.pop()
                    pairs.append([j, i])

        pairs.sort(key=lambda p: p[0])

        num_pairs = len(pairs)
        paired_fraction = (2 * num_pairs) / n if n > 0 else 0.0

        # Confidence estimate based on MFE magnitude
        confidence = min(1.0, max(0.1, abs(mfe) / (n * 0.5))) if n > 0 else 0.5

        # Ensemble diversity estimate
        diversity = _estimate_diversity(n, num_pairs) if n > 0 else 0.0

        return FoldResult(
            dot_bracket=dot_bracket,
            minimum_free_energy=round(mfe, 2),
            base_pairs=pairs,
            structural_confidence=round(confidence, 4),
            paired_fraction=round(paired_fraction, 4),
            ensemble_diversity=round(diversity, 4),
        )

    except Exception as e:
        logger.warning("ViennaRNA prediction failed: %s — falling back", e)
        return None


def fold_with_nussinov(sequence: str) -> FoldResult:
    """Predict secondary structure using the pure-Python Nussinov engine."""
    from app.services.structure_prediction import predict_structure

    result = predict_structure(sequence)
    return FoldResult(
        dot_bracket=result.dot_bracket,
        minimum_free_energy=result.minimum_free_energy,
        base_pairs=result.base_pairs,
        structural_confidence=result.structural_confidence,
        paired_fraction=result.paired_fraction,
        ensemble_diversity=result.ensemble_diversity,
    )


def predict_secondary_structure(sequence: str) -> FoldResult:
    """Predict RNA secondary structure with automatic engine selection.

    1. Try ViennaRNA (if installed)
    2. Fall back to pure-Python Nussinov
    """
    result = fold_with_vienna(sequence)
    if result is not None:
        return result
    return fold_with_nussinov(sequence)


def _estimate_diversity(n: int, num_pairs: int) -> float:
    """Estimate ensemble diversity from pairing density."""
    import math
    max_possible = n // 2
    if max_possible == 0:
        return 0.0
    ratio = num_pairs / max_possible
    diversity = 4.0 * ratio * (1.0 - ratio)
    length_factor = min(1.0, math.log(n) / math.log(500))
    diversity *= length_factor
    return diversity
