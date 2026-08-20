from __future__ import annotations

import random
from typing import Tuple

from app.schemas.analysis_result import AnalysisResult
from app.schemas.rna_sequence import RNASequenceInput
from app.services.rna_analysis import analyze_sequence

_NUCLEOTIDES = ["A", "C", "G", "U"]


def generate_rna(length: int, gc_content: float) -> str:
    """Generate a random RNA sequence with a target GC-content.

    Args:
        length:   Desired sequence length (1–10_000).
        gc_content: Target fraction of G+C nucleotides (0.0–1.0).

    Returns:
        An uppercase RNA string composed only of A, C, G, U.
    """
    if length < 1:
        raise ValueError("Sequence length must be at least 1.")
    if length > 10_000:
        raise ValueError("Sequence length must not exceed 10_000.")
    if not 0.0 <= gc_content <= 1.0:
        raise ValueError("GC-content must be between 0.0 and 1.0.")

    gc_count = round(length * gc_content)
    au_count = length - gc_count

    pool: list[str] = (
        random.choices(["G", "C"], k=gc_count) +
        random.choices(["A", "U"], k=au_count)
    )
    random.shuffle(pool)
    return "".join(pool)


def generate_and_analyze(length: int, gc_content: float) -> AnalysisResult:
    """Generate an RNA sequence and run the full analysis pipeline on it.

    The result is automatically cached in the analysis registry so that
    structure, candidate, and export pages can reference it.

    Returns a fully populated AnalysisResult with real computed metrics.
    """
    sequence = generate_rna(length, gc_content)
    name = f"DeNovo-GC{gc_content:.2f}-L{length}"

    input_data = RNASequenceInput(sequence=sequence, name=name)
    return analyze_sequence(input_data)


def validate_rna(sequence: str) -> Tuple[bool, str]:
    """Check that a string is a valid RNA sequence."""
    if not sequence:
        return False, "Sequence is empty."
    allowed = set("ACGU")
    found = set(sequence.upper())
    invalid = found - allowed
    if invalid:
        return False, f"Invalid characters: {invalid}"
    return True, ""
