from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field


class ExportRequest(BaseModel):
    """Request to export analysis results in a specified format."""

    analysis_ids: list[str] = Field(
        ...,
        min_length=1,
        description="List of analysis IDs to include in the export.",
    )
    format: Literal["pdf", "fasta", "json"] = Field(
        ...,
        description="Output format: 'pdf' for report, 'fasta' for sequence file, 'json' for raw data.",
    )
    include_visualizations: bool = Field(
        default=True,
        description="Whether to include structure visualizations (PDF only).",
    )
    include_interpretation: bool = Field(
        default=True,
        description="Whether to include genomic interpretation text.",
    )
