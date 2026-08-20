"""
Standalone RNA analysis endpoint module.

Provides the POST /analyze endpoint returning the exact response
format specified for the genomic analysis pipeline.
"""

from __future__ import annotations

from fastapi import APIRouter, HTTPException

from app.schemas.analyze_payloads import AnalyzeRequest, AnalyzeResponse
from app.services.analyze_engine import run_analysis

router = APIRouter(tags=["RNA Analysis"])


@router.post("/analyze", response_model=AnalyzeResponse)
async def analyze_candidate(payload: AnalyzeRequest):
    """
    Submit an RNA sequence for full genomic analysis.

    Performs sequence normalization, validation, secondary structure
    prediction, thermodynamic stability scoring, and generates a
    molecular interpretation. Returns metrics in a flat production format.
    """
    try:
        return run_analysis(payload)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Analysis pipeline error: {str(e)}",
        )
