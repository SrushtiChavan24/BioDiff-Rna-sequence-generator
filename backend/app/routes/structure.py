"""
RNA Structure Routes.

Handles structure prediction queries and generates SVG structure diagrams.
"""

from __future__ import annotations

from fastapi import APIRouter, HTTPException, Response
from fastapi.responses import HTMLResponse

from app.schemas.rna_sequence import RNASequenceInput
from app.schemas.analysis_result import StructurePrediction
from app.services.structure_prediction import predict_structure
from app.services.rna_analysis import get_cached_analysis
from app.visualization.structure_svg import generate_structure_svg
from app.services.structure_generation import generate_true_fold_svg
from app.utils.validators import SequenceValidationError

router = APIRouter(prefix="/structure", tags=["RNA Structure"])


@router.post("/predict", response_model=StructurePrediction)
async def predict_secondary_structure(payload: RNASequenceInput):
    """
    Predict the secondary structure (dot-bracket notation) for a sequence
    without performing a full analysis lifecycle.
    """
    try:
        struct_res = predict_structure(payload.sequence)
        return StructurePrediction(
            dot_bracket=struct_res.dot_bracket,
            minimum_free_energy=struct_res.minimum_free_energy,
            base_pairs=struct_res.base_pairs,
            structural_confidence=struct_res.structural_confidence,
            paired_fraction=struct_res.paired_fraction,
        )
    except SequenceValidationError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"An error occurred during structure prediction: {str(e)}",
        )


@router.get("/{analysis_id}/svg")
async def get_structure_svg(analysis_id: str, width: int = 600, height: int = 500):
    """
    Retrieve the true RNA secondary structure SVG for a past analysis.

    Returns a publication-quality stem-loop topology diagram
    as raw SVG XML with correct Content-Type.
    """
    cached = get_cached_analysis(analysis_id)
    if not cached:
        raise HTTPException(
            status_code=404,
            detail=f"Analysis with ID '{analysis_id}' not found.",
        )

    svg_content = generate_true_fold_svg(
        sequence=cached.sequence,
        dot_bracket=cached.structure.dot_bracket,
        base_pairs=cached.structure.base_pairs,
        width=width,
        height=height,
        title=cached.sequence_name,
    )

    return Response(content=svg_content, media_type="image/svg+xml")


@router.get("/{analysis_id}/mountain")
async def get_mountain_plot_svg(analysis_id: str, width: int = 600, height: int = 500):
    """
    Retrieve the classic mountain-plot SVG for a past analysis.
    """
    cached = get_cached_analysis(analysis_id)
    if not cached:
        raise HTTPException(
            status_code=404,
            detail=f"Analysis with ID '{analysis_id}' not found.",
        )

    svg_content = generate_structure_svg(
        sequence=cached.sequence,
        dot_bracket=cached.structure.dot_bracket,
        base_pairs=cached.structure.base_pairs,
        width=width,
        height=height,
        title=cached.sequence_name,
    )

    return Response(content=svg_content, media_type="image/svg+xml")
