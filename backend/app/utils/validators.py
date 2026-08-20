"""
RNA sequence validation utilities.

Provides validation functions used by routes and services to ensure
sequence data integrity before computation.
"""

from __future__ import annotations

import re

from app.config import settings


class SequenceValidationError(ValueError):
    """Raised when an RNA sequence fails validation."""

    pass


def validate_rna_sequence(sequence: str) -> str:
    """
    Validate and normalize an RNA sequence string.

    - Strips whitespace and newlines
    - Converts to uppercase
    - Converts T → U (DNA to RNA)
    - Rejects any character not in {A, C, G, U}

    Returns the cleaned sequence.
    Raises SequenceValidationError on invalid input.
    """
    if not sequence or not sequence.strip():
        raise SequenceValidationError("Sequence cannot be empty.")

    cleaned = sequence.upper().replace(" ", "").replace("\n", "").replace("\r", "")
    cleaned = cleaned.replace("T", "U")

    if not re.fullmatch(r"[ACGU]+", cleaned):
        invalid = sorted(set(re.findall(r"[^ACGU]", cleaned)))
        raise SequenceValidationError(
            f"Invalid nucleotide characters: {invalid}. "
            f"RNA sequences must contain only A, C, G, U."
        )

    if len(cleaned) < settings.MIN_SEQUENCE_LENGTH:
        raise SequenceValidationError(
            f"Sequence too short ({len(cleaned)} nt). "
            f"Minimum length: {settings.MIN_SEQUENCE_LENGTH} nt."
        )

    if len(cleaned) > settings.MAX_SEQUENCE_LENGTH:
        raise SequenceValidationError(
            f"Sequence too long ({len(cleaned)} nt). "
            f"Maximum length: {settings.MAX_SEQUENCE_LENGTH} nt."
        )

    return cleaned


def validate_batch_size(count: int) -> None:
    """Validate that the batch size is within limits."""
    if count > settings.MAX_BATCH_SIZE:
        raise SequenceValidationError(
            f"Batch size {count} exceeds maximum of {settings.MAX_BATCH_SIZE}."
        )
    if count < 1:
        raise SequenceValidationError("Batch must contain at least 1 sequence.")
