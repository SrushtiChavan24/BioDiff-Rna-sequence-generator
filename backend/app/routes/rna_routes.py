from __future__ import annotations

import uuid

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

from app.services.rna_generator import generate_and_analyze
from app.services.rna_analysis import get_all_cached_analyses

router = APIRouter(prefix="/rna", tags=["RNA Generation"])


class GenerateRequest(BaseModel):
    length: int = Field(default=50, ge=1, le=10_000, description="Desired sequence length")
    gc_content: float = Field(default=0.5, ge=0.0, le=1.0, description="Target GC fraction")


class GenerateResponse(BaseModel):
    candidate_id: str
    analysis_id: str
    sequence: str
    length: int
    gc_content: float
    mfe: float
    paired_ratio: float
    therapeutic_score: float
    structure_confidence: float
    dot_bracket: str
    interpretation: str
    ranking: int


def _compute_ranking(score: float) -> int:
    """Determine rank among all cached candidates by therapeutic score."""
    all_cached = get_all_cached_analyses()
    if not all_cached:
        return 1
    sorted_scores = sorted(
        [a.stability.therapeutic_stability_score for a in all_cached],
        reverse=True,
    )
    for rank, s in enumerate(sorted_scores, start=1):
        if score >= s:
            return rank
    return len(sorted_scores) + 1


@router.post("/generate", response_model=GenerateResponse)
async def rna_generate(payload: GenerateRequest):
    """Generate an RNA sequence and run the full computational genomics pipeline.

    Returns the sequence along with real computed metrics:
    - MFE via Nussinov DP with thermodynamic parameters
    - Therapeutic stability score (0–100)
    - Dot-bracket secondary structure
    - Biotech-grade genomic interpretation
    - Candidate ranking within the cohort
    """
    try:
        result = generate_and_analyze(payload.length, payload.gc_content)
        score = result.stability.therapeutic_stability_score
        return GenerateResponse(
            candidate_id=f"RNA-{uuid.uuid4().hex[:6].upper()}",
            analysis_id=result.analysis_id,
            sequence=result.sequence,
            length=result.length,
            gc_content=round(result.stability.gc_content, 4),
            mfe=round(result.structure.minimum_free_energy, 2),
            paired_ratio=round(result.structure.paired_fraction, 4),
            therapeutic_score=round(score, 1),
            structure_confidence=round(result.structure.structural_confidence * 100, 1),
            dot_bracket=result.structure.dot_bracket,
            interpretation=result.genomic_interpretation,
            ranking=_compute_ranking(score),
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
