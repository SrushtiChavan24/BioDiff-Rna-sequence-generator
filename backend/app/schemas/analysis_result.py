from __future__ import annotations

from pydantic import BaseModel, Field


class StructurePrediction(BaseModel):
    """RNA secondary structure prediction result."""

    dot_bracket: str = Field(
        ...,
        description="Dot-bracket notation of the predicted minimum free energy structure.",
    )
    minimum_free_energy: float = Field(
        ...,
        description="Minimum free energy of the predicted structure in kcal/mol.",
    )
    base_pairs: list[list[int]] = Field(
        ...,
        description="List of base pair positions [i, j] (0-indexed).",
    )
    structural_confidence: float = Field(
        ...,
        ge=0.0,
        le=1.0,
        description=(
            "Confidence in the predicted structure, derived from ensemble "
            "analysis. 1.0 indicates a single dominant fold; lower values "
            "indicate structural heterogeneity."
        ),
    )
    paired_fraction: float = Field(
        ...,
        ge=0.0,
        le=1.0,
        description="Fraction of nucleotides participating in base pairs.",
    )


class NucleotideComposition(BaseModel):
    """Nucleotide frequency breakdown."""

    A: float = Field(..., description="Adenine fraction")
    C: float = Field(..., description="Cytosine fraction")
    G: float = Field(..., description="Guanine fraction")
    U: float = Field(..., description="Uracil fraction")


class StabilityMetrics(BaseModel):
    """Thermodynamic and structural stability assessment."""

    therapeutic_stability_score: float = Field(
        ...,
        ge=0.0,
        le=100.0,
        description=(
            "Composite stability score (0-100) integrating thermodynamic "
            "stability, structural confidence, and GC content. Higher values "
            "indicate greater suitability as a therapeutic candidate."
        ),
    )
    gc_content: float = Field(
        ...,
        ge=0.0,
        le=1.0,
        description="Guanine-cytosine content ratio.",
    )
    mfe_per_nucleotide: float = Field(
        ...,
        description="Minimum free energy normalized per nucleotide (kcal/mol/nt).",
    )
    ensemble_diversity: float = Field(
        ...,
        ge=0.0,
        description=(
            "Structural ensemble diversity — lower values indicate a single "
            "dominant conformation; higher values indicate conformational flexibility."
        ),
    )


class AnalysisResult(BaseModel):
    """Complete analysis output for a single RNA sequence."""

    analysis_id: str = Field(
        ...,
        description="Unique identifier for this analysis.",
    )
    sequence_name: str = Field(
        ...,
        description="Name or auto-generated label for this sequence.",
    )
    sequence: str = Field(
        ...,
        description="The analyzed RNA sequence.",
    )
    length: int = Field(
        ...,
        description="Nucleotide count.",
    )
    composition: NucleotideComposition = Field(
        ...,
        description="Nucleotide composition breakdown.",
    )
    structure: StructurePrediction = Field(
        ...,
        description="Predicted secondary structure.",
    )
    stability: StabilityMetrics = Field(
        ...,
        description="Thermodynamic stability metrics.",
    )
    genomic_interpretation: str = Field(
        ...,
        description=(
            "Human-readable, medically oriented interpretation of the "
            "analysis results."
        ),
    )
    timestamp: str = Field(
        ...,
        description="ISO 8601 timestamp of analysis completion.",
    )


class BatchAnalysisResult(BaseModel):
    """Result container for batch analysis."""

    results: list[AnalysisResult]
    total_sequences: int
    successful: int
    failed: int


class SimilarityResult(BaseModel):
    """Molecular similarity landscape computation result."""

    similarity_matrix: list[list[float]] = Field(
        ...,
        description="N×N similarity matrix (values 0–1, 1 = identical).",
    )
    labels: list[str] = Field(
        ...,
        description="Sequence labels corresponding to matrix rows/columns.",
    )
    analysis_ids: list[str] = Field(
        ...,
        description="Analysis IDs corresponding to matrix rows/columns.",
    )


class StructuralFamily(BaseModel):
    """A detected structural family (cluster) of RNA sequences."""

    family_id: int
    family_label: str
    member_ids: list[str]
    centroid_id: str
    avg_similarity: float


class ProjectionPoint(BaseModel):
    """A single point in the structural relationship projection."""

    analysis_id: str
    label: str
    x: float
    y: float
    family_id: int


class ClusteringResult(BaseModel):
    """Structural family detection and projection result."""

    families: list[StructuralFamily]
    projection: list[ProjectionPoint] = Field(
        ...,
        description="2D structural relationship projection coordinates.",
    )
    total_families: int
