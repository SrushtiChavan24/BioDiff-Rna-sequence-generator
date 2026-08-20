from __future__ import annotations

from pydantic import BaseModel, Field, field_validator
import re


class AnalyzeRequest(BaseModel):
    name: str = Field(
        default="unnamed_candidate",
        max_length=256,
        description="Human-readable candidate identifier.",
        examples=["Spike_vaccine_v2"],
    )
    sequence: str = Field(
        ...,
        min_length=10,
        max_length=10_000,
        description="Raw RNA nucleotide sequence (A, C, G, U / T).",
        examples=["GGGCGCGCGCCGAUGCUAGCUAGCUAGC"],
    )

    @field_validator("sequence")
    @classmethod
    def normalize_sequence(cls, v: str) -> str:
        cleaned = v.upper().replace(" ", "").replace("\n", "").replace("\r", "")
        cleaned = cleaned.replace("T", "U")
        if not re.fullmatch(r"[ACGU]+", cleaned):
            invalid = sorted(set(re.findall(r"[^ACGU]", cleaned)))
            raise ValueError(
                f"Invalid nucleotide characters detected: {invalid}. "
                f"Only A, C, G, U (or T) are permitted."
            )
        return cleaned


class AnalyzeResponse(BaseModel):
    candidate_id: str = Field(..., description="Unique candidate identifier.", examples=["RNA-001"])
    analysis_id: str = Field(..., description="Analysis cache ID for structure page navigation.", examples=["ANL-A1B2C3D4"])
    name: str = Field(..., description="Candidate name.", examples=["Spike_vaccine_v2"])
    sequence: str = Field(..., description="Normalized RNA sequence.")
    length: int = Field(..., description="Nucleotide count.")
    gc_content: float = Field(..., description="GC content percentage (0–100).", examples=[63.3])
    au_content: float = Field(..., description="AU content percentage (0–100).", examples=[36.7])
    mfe: float = Field(..., description="Minimum free energy in kcal/mol.", examples=[-18.2])
    paired_ratio: float = Field(..., description="Fraction of nucleotides in base pairs (0–1).", examples=[0.71])
    therapeutic_score: float = Field(..., description="Composite therapeutic stability score (0–100).", examples=[84.5])
    structure_confidence: float = Field(..., description="Confidence in predicted fold (0–100).", examples=[91.2])
    dot_bracket: str = Field(..., description="Dot-bracket secondary structure notation.")
    base_pairs: list[list[int]] = Field(..., description="List of base pair positions [i, j] (0-indexed).")
    interpretation: str = Field(..., description="Molecular and therapeutic interpretation.")
