"""
RNA Analysis Orchestrator Service.

Coordinates structure prediction, stability metrics, composition analysis,
and genomic interpretation to construct a unified AnalysisResult.
Manages an in-memory cache of results.
"""

from __future__ import annotations

from datetime import datetime
import uuid

from app.schemas.rna_sequence import RNASequenceInput
from app.schemas.analysis_result import (
    AnalysisResult,
    NucleotideComposition,
    StructurePrediction,
    StabilityMetrics,
)
from app.utils.sequence_utils import compute_nucleotide_frequencies, compute_gc_content
from app.services.structure_prediction import predict_structure
from app.services.stability import compute_stability_metrics
from app.interpretations.genomic_interpreter import interpret_analysis

# In-memory database of analysis results
# Keys are analysis_id, values are AnalysisResult instances
_ANALYSIS_CACHE: dict[str, AnalysisResult] = {}


def get_cached_analysis(analysis_id: str) -> AnalysisResult | None:
    """Retrieve an analysis result from the cache by ID."""
    return _ANALYSIS_CACHE.get(analysis_id)


def get_all_cached_analyses() -> list[AnalysisResult]:
    """Retrieve all cached analysis results, sorted by timestamp descending."""
    return sorted(
        _ANALYSIS_CACHE.values(),
        key=lambda x: x.timestamp,
        reverse=True,
    )


def clear_cached_analyses() -> None:
    """Clear all results from the cache."""
    _ANALYSIS_CACHE.clear()


def analyze_sequence(sequence_input: RNASequenceInput) -> AnalysisResult:
    """
    Orchestrate full RNA sequence analysis.

    Runs composition analysis, structure prediction, stability scoring,
    and clinical interpretation. Stores the result in the cache.

    Args:
        sequence_input: Validated input containing the RNA sequence.

    Returns:
        Complete AnalysisResult.
    """
    sequence = sequence_input.sequence
    name = sequence_input.name or f"RNA-Seq-{uuid.uuid4().hex[:6].upper()}"

    # 1. Compute Composition
    frequencies = compute_nucleotide_frequencies(sequence)
    composition = NucleotideComposition(
        A=frequencies["A"],
        C=frequencies["C"],
        G=frequencies["G"],
        U=frequencies["U"],
    )

    # 2. Predict Secondary Structure
    struct_result = predict_structure(sequence)
    structure = StructurePrediction(
        dot_bracket=struct_result.dot_bracket,
        minimum_free_energy=struct_result.minimum_free_energy,
        base_pairs=struct_result.base_pairs,
        structural_confidence=struct_result.structural_confidence,
        paired_fraction=struct_result.paired_fraction,
    )

    # 3. Compute Stability Metrics
    stability_result = compute_stability_metrics(
        sequence=sequence,
        mfe=struct_result.minimum_free_energy,
        structural_confidence=struct_result.structural_confidence,
        paired_fraction=struct_result.paired_fraction,
        ensemble_diversity=struct_result.ensemble_diversity,
    )
    stability = StabilityMetrics(
        therapeutic_stability_score=stability_result.therapeutic_stability_score,
        gc_content=stability_result.gc_content,
        mfe_per_nucleotide=stability_result.mfe_per_nucleotide,
        ensemble_diversity=stability_result.ensemble_diversity,
    )

    # 4. Generate Genomic Interpretation
    interpretation = interpret_analysis(
        sequence_name=name,
        sequence_length=len(sequence),
        gc_content=stability_result.gc_content,
        mfe=struct_result.minimum_free_energy,
        mfe_per_nt=stability_result.mfe_per_nucleotide,
        structural_confidence=struct_result.structural_confidence,
        paired_fraction=struct_result.paired_fraction,
        ensemble_diversity=stability_result.ensemble_diversity,
        therapeutic_stability_score=stability_result.therapeutic_stability_score,
        num_base_pairs=len(struct_result.base_pairs),
    )

    # 5. Construct full result
    analysis_id = f"ANL-{uuid.uuid4().hex[:8].upper()}"
    result = AnalysisResult(
        analysis_id=analysis_id,
        sequence_name=name,
        sequence=sequence,
        length=len(sequence),
        composition=composition,
        structure=structure,
        stability=stability,
        genomic_interpretation=interpretation,
        timestamp=datetime.utcnow().isoformat() + "Z",
    )

    # Cache the result
    _ANALYSIS_CACHE[analysis_id] = result
    return result


def analyze_batch(batch_input: list[RNASequenceInput]) -> list[AnalysisResult]:
    """Analyze a batch of RNA sequences and return the results."""
    results = []
    for seq_in in batch_input:
        try:
            results.append(analyze_sequence(seq_in))
        except Exception as e:
            # In a production system, we'd log this and potentially include a failure in the response.
            # Here we let it propagate or skip. Let's let it propagate for safety unless handling is preferred.
            raise e
    return results
