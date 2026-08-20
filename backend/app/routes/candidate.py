from __future__ import annotations

from fastapi import APIRouter, HTTPException, Query

from app.schemas.analysis_result import (
    AnalysisResult,
    ClusteringResult,
    ProjectionPoint,
    SimilarityResult,
    StructuralFamily,
)
from app.schemas.response import BaseResponse
from app.services.clustering import detect_structural_families
from app.services.rna_analysis import get_all_cached_analyses, get_cached_analysis
from app.services.similarity import compute_similarity_matrix

router = APIRouter(prefix="/candidates", tags=["Therapeutic Candidates"])
legacy_router = APIRouter(tags=["Therapeutic Candidates (Legacy)"])


@router.get("", response_model=BaseResponse[list[AnalysisResult]])
@router.get("/registry", response_model=BaseResponse[list[AnalysisResult]])
@legacy_router.get("/candidates", response_model=BaseResponse[list[AnalysisResult]], include_in_schema=False)
async def list_candidate_registry(
    sort_by: str = Query(
        default="stability",
        description="Sort field: 'stability', 'length', 'gc', or 'name'.",
    ),
    ascending: bool = Query(default=False, description="Sort ascending if true."),
):
    """
    Retrieve all profiled therapeutic candidates from the registry.

    Results are sorted according to the specified field and order, enabling
    rapid identification of high-priority candidates for further analysis.
    """
    analyses = get_all_cached_analyses()

    sort_key_map = {
        "stability": lambda a: a.stability.therapeutic_stability_score,
        "length": lambda a: a.length,
        "gc": lambda a: a.stability.gc_content,
        "name": lambda a: a.sequence_name.lower(),
    }
    key_fn = sort_key_map.get(sort_by, sort_key_map["stability"])
    analyses.sort(key=key_fn, reverse=not ascending)

    return BaseResponse(
        success=True,
        message=f"Found {len(analyses)} candidate(s) in registry.",
        data=analyses,
    )


@router.get("/rankings", response_model=BaseResponse[list[AnalysisResult]])
async def get_candidate_rankings(
    min_score: float = Query(
        default=0.0,
        ge=0.0,
        le=100.0,
        description="Minimum Therapeutic Stability Score threshold.",
    ),
    max_results: int = Query(
        default=50, ge=1, le=100, description="Maximum number of candidates to return."
    ),
):
    """
    Rank profiled candidates by Therapeutic Stability Score.

    Filters to candidates meeting the minimum score threshold and returns
    them sorted from highest to lowest stability index.  Provides a rapid
    triage view for therapeutic candidate selection.
    """
    analyses = get_all_cached_analyses()
    filtered = [
        a
        for a in analyses
        if a.stability.therapeutic_stability_score >= min_score
    ]
    filtered.sort(
        key=lambda a: a.stability.therapeutic_stability_score, reverse=True
    )
    ranked = filtered[:max_results]

    return BaseResponse(
        success=True,
        message=(
            f"Ranked {len(ranked)} candidate(s) meeting ≥ {min_score}/100 "
            f"threshold (of {len(analyses)} total)."
        ),
        data=ranked,
    )


@router.get(
    "/{analysis_id}", response_model=BaseResponse[AnalysisResult]
)
@legacy_router.get(
    "/candidate/{analysis_id}", response_model=BaseResponse[AnalysisResult], include_in_schema=False
)
async def get_candidate_detail(analysis_id: str):
    """
    Retrieve full analysis details for a specific therapeutic candidate.
    """
    result = get_cached_analysis(analysis_id)
    if not result:
        raise HTTPException(
            status_code=404,
            detail=f"Candidate with ID '{analysis_id}' not found in registry.",
        )
    return BaseResponse(
        success=True,
        message=f"Candidate '{result.sequence_name}' retrieved.",
        data=result,
    )


@router.post("/compare", response_model=BaseResponse[SimilarityResult])
async def compare_candidates(analysis_ids: list[str]):
    """
    Compute the molecular similarity landscape across selected candidates.

    Produces an N×N similarity matrix combining sequence, structural, and
    compositional metrics.
    """
    if len(analysis_ids) < 2:
        raise HTTPException(
            status_code=400,
            detail="At least two candidates must be selected for comparison.",
        )

    analyses: list[AnalysisResult] = []
    for aid in analysis_ids:
        cached = get_cached_analysis(aid)
        if cached:
            analyses.append(cached)

    if len(analyses) < 2:
        raise HTTPException(
            status_code=400,
            detail="Fewer than two valid analysis records found for the provided IDs.",
        )

    sequences = [a.sequence for a in analyses]
    structures = [a.structure.dot_bracket for a in analyses]
    labels = [a.sequence_name for a in analyses]

    matrix = compute_similarity_matrix(sequences, structures)

    return BaseResponse(
        success=True,
        message=f"Similarity landscape computed for {len(analyses)} candidates.",
        data=SimilarityResult(
            similarity_matrix=matrix.tolist(),
            labels=labels,
            analysis_ids=[a.analysis_id for a in analyses],
        ),
    )


@router.post("/cluster", response_model=BaseResponse[ClusteringResult])
async def detect_candidate_families(
    analysis_ids: list[str],
    n_families: int | None = Query(
        default=None,
        description="Number of structural families to detect (auto-determined if omitted).",
    ),
):
    """
    Detect candidate structural families among the selected candidates.

    Projects sequences into a 2D "Sequence Relationship Space" and assigns
    each to a structural family based on multi-metric similarity.
    """
    if len(analysis_ids) < 2:
        raise HTTPException(
            status_code=400,
            detail="At least two candidates are required for family detection.",
        )

    analyses: list[AnalysisResult] = []
    for aid in analysis_ids:
        cached = get_cached_analysis(aid)
        if cached:
            analyses.append(cached)

    if len(analyses) < 2:
        raise HTTPException(
            status_code=400,
            detail="Fewer than two valid analysis records found.",
        )

    sequences = [a.sequence for a in analyses]
    structures = [a.structure.dot_bracket for a in analyses]
    labels = [a.sequence_name for a in analyses]

    similarity_matrix = compute_similarity_matrix(sequences, structures)
    clustering_out = detect_structural_families(
        similarity_matrix=similarity_matrix,
        labels=labels,
        n_families=n_families,
    )

    families_models = []
    for fam in clustering_out.families:
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
    for pt in clustering_out.projection:
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

    return BaseResponse(
        success=True,
        message=f"Detected {len(families_models)} structural family(ies) across {len(analyses)} candidates.",
        data=ClusteringResult(
            families=families_models,
            projection=projection_models,
            total_families=len(families_models),
        ),
    )
