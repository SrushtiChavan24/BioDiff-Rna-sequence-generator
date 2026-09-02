"""
Insights API Routes.

Provides endpoints for computational genomics insights — top candidate
identification, comparative cohort summaries, and scientific interpretation
— without modifying existing analysis, structure, or export endpoints.
"""

from __future__ import annotations

from fastapi import APIRouter, HTTPException

from app.interpretations.comparative_analysis import generate_comparative_summary
from app.interpretations.scientific_interpreter import generate_scientific_interpretation
from app.services.rna_analysis import get_all_cached_analyses, get_cached_analysis

router = APIRouter(prefix="/insights", tags=["Genomics Insights"])


@router.get("/top-candidate")
async def get_top_candidate_insight():
    """
    Retrieve the top-ranked therapeutic candidate with scientific interpretation.

    Returns the highest-scoring candidate by Therapeutic Stability Index
    along with a professional biotech-grade scientific interpretation.
    """
    analyses = get_all_cached_analyses()
    if not analyses:
        raise HTTPException(
            status_code=404,
            detail="No analyzed candidates available. Submit a sequence for analysis first.",
        )

    top = max(
        analyses,
        key=lambda a: a.stability.therapeutic_stability_score,
    )

    interpretation = generate_scientific_interpretation(top)

    return {
        "success": True,
        "data": {
            "analysis_id": top.analysis_id,
            "sequence_name": top.sequence_name,
            "therapeutic_score": top.stability.therapeutic_stability_score,
            "gc_content": top.stability.gc_content,
            "mfe": top.structure.minimum_free_energy,
            "structural_confidence": top.structure.structural_confidence,
            "paired_fraction": top.structure.paired_fraction,
            "ensemble_diversity": top.stability.ensemble_diversity,
            "interpretation": interpretation,
        },
    }


@router.get("/comparative-summary")
async def get_comparative_summary():
    """
    Generate a structured comparative analysis across all profiled candidates.

    Returns cohort-level metrics, top-performer identification, strengths,
    weaknesses, and a professional recommendation.
    """
    analyses = get_all_cached_analyses()
    if not analyses:
        raise HTTPException(
            status_code=404,
            detail="No analyzed candidates available for comparison.",
        )

    if len(analyses) < 2:
        # Single candidate — return individual insight instead
        single = analyses[0]
        interpretation = generate_scientific_interpretation(single)
        return {
            "success": True,
            "data": {
                "cohort_size": 1,
                "top_candidate": single.sequence_name,
                "mean_score": single.stability.therapeutic_stability_score,
                "summary": f"Single candidate '{single.sequence_name}' analyzed. "
                           f"Submit additional sequences for comparative cohort analysis.",
                "strengths": [],
                "weaknesses": [],
                "recommendation": "Profile additional candidates to enable comparative assessment.",
                "individual_interpretation": interpretation,
            },
        }

    summary = generate_comparative_summary(analyses)
    return {
        "success": True,
        "data": summary,
    }


@router.get("/interpret/{analysis_id}")
async def get_candidate_interpretation(analysis_id: str):
    """
    Generate a detailed scientific interpretation for a specific candidate.
    """
    cached = get_cached_analysis(analysis_id)
    if not cached:
        raise HTTPException(
            status_code=404,
            detail=f"Analysis with ID '{analysis_id}' not found.",
        )

    interpretation = generate_scientific_interpretation(cached)

    return {
        "success": True,
        "data": {
            "analysis_id": cached.analysis_id,
            "sequence_name": cached.sequence_name,
            "interpretation": interpretation,
        },
    }
