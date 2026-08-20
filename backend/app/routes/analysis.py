"""
RNA Analysis Route Module.

Handles sequence submission, batch analysis, history retrieval,
molecular similarity landscape, and structural family detection.
"""

from __future__ import annotations

from fastapi import APIRouter, HTTPException, Query

from app.schemas.rna_sequence import RNASequenceInput, RNABatchInput
from app.schemas.analysis_result import (
    AnalysisResult,
    BatchAnalysisResult,
    SimilarityResult,
    ClusteringResult,
    StructuralFamily,
    ProjectionPoint,
)
from app.services.rna_analysis import (
    analyze_sequence,
    analyze_batch,
    get_cached_analysis,
    get_all_cached_analyses,
    clear_cached_analyses,
)
from app.services.similarity import compute_similarity_matrix
from app.services.clustering import detect_structural_families
from app.utils.validators import validate_batch_size, SequenceValidationError

router = APIRouter(prefix="/analysis", tags=["RNA Analysis"])


@router.post("/sequence", response_model=AnalysisResult)
async def analyze_single_sequence(payload: RNASequenceInput):
    """
    Submit a single RNA sequence for deep analysis.

    Performs structural prediction, thermodynamic stability evaluation,
    and produces a biotech-grade clinical interpretation.
    """
    try:
        result = analyze_sequence(payload)
        return result
    except SequenceValidationError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"An error occurred during sequence analysis: {str(e)}",
        )


@router.post("/batch", response_model=BatchAnalysisResult)
async def analyze_sequence_batch(payload: RNABatchInput):
    """
    Submit a batch of up to 50 RNA sequences for analysis.

    Returns the analysis results for each successfully processed sequence.
    """
    try:
        validate_batch_size(len(payload.sequences))
        results = analyze_batch(payload.sequences)
        return BatchAnalysisResult(
            results=results,
            total_sequences=len(payload.sequences),
            successful=len(results),
            failed=0,
        )
    except SequenceValidationError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"An error occurred during batch analysis: {str(e)}",
        )


@router.get("/history", response_model=list[AnalysisResult])
async def get_analysis_history():
    """Retrieve all cached analysis records sorted by recency."""
    return get_all_cached_analyses()


@router.delete("/history", tags=["System"])
async def clear_analysis_history():
    """Clear all analysis records from the cache."""
    clear_cached_analyses()
    return {"message": "Analysis history successfully cleared."}


@router.get("/{analysis_id}", response_model=AnalysisResult)
async def get_analysis_details(analysis_id: str):
    """Retrieve details of a specific past analysis by ID."""
    result = get_cached_analysis(analysis_id)
    if not result:
        raise HTTPException(
            status_code=404,
            detail=f"Analysis with ID '{analysis_id}' not found.",
        )
    return result


@router.post("/similarity", response_model=SimilarityResult)
async def compute_similarity_landscape(analysis_ids: list[str]):
    """
    Compute the molecular similarity landscape across selected analyses.

    Produces an N×N matrix representing sequence, structure, and
    composition similarities — the molecular similarity landscape.
    """
    if not analysis_ids:
        raise HTTPException(
            status_code=400,
            detail="At least one analysis ID must be specified.",
        )

    analyses: list[AnalysisResult] = []
    missing_ids = []
    for aid in analysis_ids:
        cached = get_cached_analysis(aid)
        if cached:
            analyses.append(cached)
        else:
            missing_ids.append(aid)

    if missing_ids:
        raise HTTPException(
            status_code=404,
            detail=f"The following analysis records could not be found: {missing_ids}",
        )

    if len(analyses) == 0:
        raise HTTPException(
            status_code=400,
            detail="No valid analysis records selected.",
        )

    sequences = [a.sequence for a in analyses]
    structures = [a.structure.dot_bracket for a in analyses]
    labels = [a.sequence_name for a in analyses]

    matrix = compute_similarity_matrix(sequences, structures)

    return SimilarityResult(
        similarity_matrix=matrix.tolist(),
        labels=labels,
        analysis_ids=[a.analysis_id for a in analyses],
    )


@router.post("/families", response_model=ClusteringResult)
@router.post("/clustering", response_model=ClusteringResult, include_in_schema=False)
async def detect_structural_families_endpoint(
    analysis_ids: list[str],
    n_families: int | None = Query(
        default=None,
        description="Number of structural families to identify (auto-determined if omitted).",
    ),
):
    """
    Detect structural families among selected sequences and project them into
    a conformational relationship space.

    Uses multi-metric molecular similarity to group sequences into fold-based
    structural families and generate a 2D therapeutic topology mapping.
    """
    if not analysis_ids:
        raise HTTPException(
            status_code=400,
            detail="At least one analysis ID must be specified.",
        )

    analyses: list[AnalysisResult] = []
    for aid in analysis_ids:
        cached = get_cached_analysis(aid)
        if cached:
            analyses.append(cached)

    if len(analyses) == 0:
        raise HTTPException(
            status_code=400,
            detail="No valid analysis records selected.",
        )

    sequences = [a.sequence for a in analyses]
    structures = [a.structure.dot_bracket for a in analyses]
    labels = [a.sequence_name for a in analyses]

    similarity_matrix = compute_similarity_matrix(sequences, structures)

    family_out = detect_structural_families(
        similarity_matrix=similarity_matrix,
        labels=labels,
        n_families=n_families,
    )

    family_labels = {
        fam.family_id: fam.family_label for fam in family_out.families
    }

    families_models = []
    for fam in family_out.families:
        member_ids = [analyses[idx].analysis_id for idx in fam.member_indices]
        centroid_id = analyses[fam.centroid_index].analysis_id
        families_models.append(
            StructuralFamily(
                family_id=fam.family_id,
                family_label=fam.family_label,
                member_ids=member_ids,
                centroid_id=centroid_id,
                avg_similarity=fam.avg_similarity,
            )
        )

    projection_models = []
    for pt in family_out.projection:
        idx = pt.index
        projection_models.append(
            ProjectionPoint(
                analysis_id=analyses[idx].analysis_id,
                label=analyses[idx].sequence_name,
                x=pt.x,
                y=pt.y,
                family_id=pt.family_id,
            )
        )

    return ClusteringResult(
        families=families_models,
        projection=projection_models,
        total_families=len(families_models),
    )
