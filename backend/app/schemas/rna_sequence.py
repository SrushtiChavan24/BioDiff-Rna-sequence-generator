from __future__ import annotations

import re

from pydantic import BaseModel, Field, field_validator


class RNASequenceInput(BaseModel):
    """Single RNA sequence submission for analysis."""

    sequence: str = Field(
        ...,
        min_length=10,
        max_length=10_000,
        description="Raw RNA nucleotide sequence containing only A, C, G, U characters.",
        examples=["AUGCGAUUCGAUAGCUAGCUAGCUAGCUAGCU"],
    )
    name: str | None = Field(
        default=None,
        max_length=256,
        description="Optional human-readable identifier for this sequence.",
        examples=["tRNA-Ala-Human"],
    )
    description: str | None = Field(
        default=None,
        max_length=1024,
        description="Optional description or notes about the sequence origin.",
    )

    @field_validator("sequence")
    @classmethod
    def validate_rna_characters(cls, v: str) -> str:
        """Normalize and validate the RNA sequence."""
        v = v.upper().replace(" ", "").replace("\n", "").replace("\r", "")
        v = v.replace("T", "U")
        if not re.fullmatch(r"[ACGU]+", v):
            invalid_chars = set(re.findall(r"[^ACGU]", v))
            raise ValueError(
                f"Invalid nucleotide characters detected: {invalid_chars}. "
                f"RNA sequences must contain only A, C, G, U."
            )
        return v


class RNABatchInput(BaseModel):
    """Batch submission of multiple RNA sequences."""

    sequences: list[RNASequenceInput] = Field(
        ...,
        min_length=1,
        max_length=50,
        description="List of RNA sequences to analyze (max 50 per batch).",
    )
