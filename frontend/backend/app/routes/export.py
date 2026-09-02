"""
Export Routes.

Provides endpoints to download analysis results as PDF reports,
FASTA files, or raw JSON data.
"""

from __future__ import annotations

import json
from fastapi import APIRouter, HTTPException, Response
from fastapi.responses import StreamingResponse
import io

from app.schemas.export_request import ExportRequest
from app.services.rna_analysis import get_cached_analysis
from app.services.report_generator import generate_pdf_report
from app.utils.sequence_utils import sequence_to_fasta

router = APIRouter(prefix="/export", tags=["Data Export"])
legacy_router = APIRouter(tags=["Data Export (Legacy)"])


@router.post("/document")
async def export_analysis_results(payload: ExportRequest):
    """
    Export selected RNA analysis records in the requested format.

    - pdf: Returns a professionally typeset PDF dossier report.
    - fasta: Returns a standard nucleotide FASTA file.
    - json: Returns a structured JSON database archive.
    """
    if not payload.analysis_ids:
        raise HTTPException(
            status_code=400,
            detail="At least one analysis ID must be provided for export.",
        )

    # Resolve analyses from cache
    analyses = []
    for aid in payload.analysis_ids:
        cached = get_cached_analysis(aid)
        if cached:
            analyses.append(cached)

    if not analyses:
        raise HTTPException(
            status_code=404,
            detail="No valid analysis records were found matching the provided IDs.",
        )

    if payload.format == "pdf":
        pdf_bytes = generate_pdf_report(analyses)
        filename = (
            f"rna_stability_dossier_{analyses[0].analysis_id}.pdf"
            if len(analyses) == 1
            else "rna_genomics_dossier.pdf"
        )
        return StreamingResponse(
            io.BytesIO(pdf_bytes),
            media_type="application/pdf",
            headers={"Content-Disposition": f"attachment; filename={filename}"},
        )

    elif payload.format == "fasta":
        fasta_content = ""
        for a in analyses:
            desc = f"length={a.length} stability={a.stability.therapeutic_stability_score}"
            fasta_content += sequence_to_fasta(
                sequence=a.sequence,
                name=a.sequence_name,
                description=desc,
            )

        filename = (
            f"rna_sequence_{analyses[0].analysis_id}.fasta"
            if len(analyses) == 1
            else "rna_sequences.fasta"
        )
        return Response(
            content=fasta_content,
            media_type="text/plain",
            headers={"Content-Disposition": f"attachment; filename={filename}"},
        )

    elif payload.format == "json":
        data = [a.model_dump() for a in analyses]
        json_content = json.dumps(data, indent=2)

        filename = (
            f"rna_analysis_{analyses[0].analysis_id}.json"
            if len(analyses) == 1
            else "rna_analyses_archive.json"
        )
        return Response(
            content=json_content,
            media_type="application/json",
            headers={"Content-Disposition": f"attachment; filename={filename}"},
        )

    else:
        raise HTTPException(
            status_code=400,
            detail=f"Unsupported export format: {payload.format}",
        )


# ---------------------------------------------------------------------------
# Format-specific export aliases
# ---------------------------------------------------------------------------


@legacy_router.post("/export/pdf", include_in_schema=False)
async def export_pdf_report(analysis_ids: list[str]):
    """Export a PDF dossier for the specified candidates."""
    req = ExportRequest(analysis_ids=analysis_ids, format="pdf")
    return await export_analysis_results(req)


@legacy_router.post("/export/json", include_in_schema=False)
async def export_json_report(analysis_ids: list[str]):
    """Export analysis data as a structured JSON archive."""
    req = ExportRequest(analysis_ids=analysis_ids, format="json")
    return await export_analysis_results(req)


@legacy_router.post("/export/fasta", include_in_schema=False)
async def export_fasta_report(analysis_ids: list[str]):
    """Export sequences in FASTA nucleotide format."""
    req = ExportRequest(analysis_ids=analysis_ids, format="fasta")
    return await export_analysis_results(req)
