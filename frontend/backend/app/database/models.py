from __future__ import annotations

import json

from sqlalchemy import Column, Integer, String, Text, Float, DateTime, func
from sqlalchemy.orm import Mapped, mapped_column

from app.database.connection import Base


class AnalysisModel(Base):
    """SQLAlchemy ORM model for persisting analysis results."""

    __tablename__ = "analyses"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    analysis_id: Mapped[str] = mapped_column(
        String(32), unique=True, nullable=False, index=True
    )
    sequence_name: Mapped[str] = mapped_column(String(256), nullable=False)
    sequence: Mapped[str] = mapped_column(Text, nullable=False)
    length: Mapped[int] = mapped_column(Integer, nullable=False)

    # Composition stored as JSON
    composition_json: Mapped[str] = mapped_column("composition", Text, nullable=False)

    # Structure fields
    dot_bracket: Mapped[str] = mapped_column(Text, nullable=False)
    minimum_free_energy: Mapped[float] = mapped_column(Float, nullable=False)
    base_pairs_json: Mapped[str] = mapped_column("base_pairs", Text, nullable=False)
    structural_confidence: Mapped[float] = mapped_column(Float, nullable=False)
    paired_fraction: Mapped[float] = mapped_column(Float, nullable=False)

    # Stability fields
    therapeutic_stability_score: Mapped[float] = mapped_column(Float, nullable=False)
    gc_content: Mapped[float] = mapped_column(Float, nullable=False)
    mfe_per_nucleotide: Mapped[float] = mapped_column(Float, nullable=False)
    ensemble_diversity: Mapped[float] = mapped_column(Float, nullable=False)

    # Interpretation
    genomic_interpretation: Mapped[str] = mapped_column(Text, nullable=False)

    # Metadata
    timestamp: Mapped[str] = mapped_column(String(32), nullable=False)
    created_at: Mapped[DateTime] = mapped_column(
        DateTime, server_default=func.now(), nullable=False
    )

    def to_dict(self) -> dict:
        """Reconstruct the AnalysisResult-compatible dictionary."""
        return {
            "analysis_id": self.analysis_id,
            "sequence_name": self.sequence_name,
            "sequence": self.sequence,
            "length": self.length,
            "composition": json.loads(self.composition_json),
            "structure": {
                "dot_bracket": self.dot_bracket,
                "minimum_free_energy": self.minimum_free_energy,
                "base_pairs": json.loads(self.base_pairs_json),
                "structural_confidence": self.structural_confidence,
                "paired_fraction": self.paired_fraction,
            },
            "stability": {
                "therapeutic_stability_score": self.therapeutic_stability_score,
                "gc_content": self.gc_content,
                "mfe_per_nucleotide": self.mfe_per_nucleotide,
                "ensemble_diversity": self.ensemble_diversity,
            },
            "genomic_interpretation": self.genomic_interpretation,
            "timestamp": self.timestamp,
        }
