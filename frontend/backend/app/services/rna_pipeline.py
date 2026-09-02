"""
Centralized RNA Computational Genomics Pipeline.

Wraps all existing backend services into a single, deterministic pipeline
that mirrors the Colab research workflow.  Every function is typed,
reusable, and independently testable.
"""

from __future__ import annotations

import uuid
from datetime import datetime
from typing import Any

import numpy as np

from app.schemas.rna_sequence import RNASequenceInput
from app.schemas.analysis_result import (
    AnalysisResult,
    NucleotideComposition,
    StabilityMetrics,
    StructurePrediction,
)
from app.services.structure_prediction import predict_structure
from app.services.stability import compute_stability_metrics
from app.services.similarity import compute_similarity_matrix as _sim_matrix
from app.services.clustering import detect_structural_families
from app.interpretations.genomic_interpreter import interpret_analysis
from app.utils.sequence_utils import (
    compute_gc_content,
    compute_nucleotide_frequencies,
    sequence_to_fasta,
)

# ---------------------------------------------------------------------------
# 1.  RNA secondary structure prediction
# ---------------------------------------------------------------------------


def generate_dot_bracket(sequence: str) -> str:
    """Return dot-bracket notation for an RNA sequence."""
    return predict_structure(sequence).dot_bracket


def compute_mfe(sequence: str) -> float:
    """Return minimum free energy (kcal/mol) for an RNA sequence."""
    return predict_structure(sequence).minimum_free_energy


def compute_paired_ratio(sequence: str) -> float:
    """Return fraction of nucleotides in base pairs."""
    return predict_structure(sequence).paired_fraction


# ---------------------------------------------------------------------------
# 2.  Therapeutic stability scoring
# ---------------------------------------------------------------------------


def compute_therapeutic_score(sequence: str) -> float:
    """Return the composite Therapeutic Stability Score (0–100)."""
    struct = predict_structure(sequence)
    metrics = compute_stability_metrics(
        sequence=sequence,
        mfe=struct.minimum_free_energy,
        structural_confidence=struct.structural_confidence,
        paired_fraction=struct.paired_fraction,
        ensemble_diversity=struct.ensemble_diversity,
    )
    return metrics.therapeutic_stability_score


# ---------------------------------------------------------------------------
# 3.  Candidate ranking
# ---------------------------------------------------------------------------


def rank_candidates(
    candidates: list[AnalysisResult],
    top_n: int | None = None,
) -> list[AnalysisResult]:
    """Rank candidates by Therapeutic Stability Score descending."""
    ranked = sorted(
        candidates,
        key=lambda c: c.stability.therapeutic_stability_score,
        reverse=True,
    )
    return ranked[:top_n] if top_n else ranked


# ---------------------------------------------------------------------------
# 4.  Similarity matrix
# ---------------------------------------------------------------------------


def compute_similarity_matrix(
    candidates: list[AnalysisResult],
    weights: dict[str, float] | None = None,
) -> np.ndarray:
    """Compute molecular similarity landscape across analysed candidates."""
    sequences = [c.sequence for c in candidates]
    structures = [c.structure.dot_bracket for c in candidates]
    return _sim_matrix(sequences, structures, weights=weights)


# ---------------------------------------------------------------------------
# 5.  PCA structural projection
# ---------------------------------------------------------------------------


def generate_pca_projection(
    candidates: list[AnalysisResult],
    n_families: int | None = None,
) -> dict[str, Any]:
    """Generate a 2D structural relationship projection with families."""
    from app.schemas.analysis_result import ClusteringResult

    sequences = [c.sequence for c in candidates]
    structures = [c.structure.dot_bracket for c in candidates]
    labels = [c.sequence_name for c in candidates]
    ids = [c.analysis_id for c in candidates]

    matrix = _sim_matrix(sequences, structures)
    out = detect_structural_families(
        similarity_matrix=matrix,
        labels=labels,
        n_families=n_families,
    )

    families = [
        {
            "family_id": f.family_id,
            "family_label": f.family_label,
            "member_ids": [ids[i] for i in f.member_indices],
            "centroid_id": ids[f.centroid_index],
            "avg_similarity": f.avg_similarity,
        }
        for f in out.families
    ]

    projection = [
        {
            "analysis_id": ids[p.index],
            "label": labels[p.index],
            "x": p.x,
            "y": p.y,
            "family_id": p.family_id,
        }
        for p in out.projection
    ]

    return {
        "families": families,
        "projection": projection,
        "total_families": out.total_families,
    }


# ---------------------------------------------------------------------------
# 6.  Full single-sequence pipeline
# ---------------------------------------------------------------------------


def run_full_analysis(
    sequence: str,
    name: str | None = None,
) -> AnalysisResult:
    """Run the complete RNA analysis pipeline on a single sequence.

    Steps:
        1. Nucleotide composition
        2. Secondary structure prediction (Nussinov DP)
        3. MFE estimation
        4. Therapeutic stability scoring
        5. Genomic interpretation

    Returns a fully populated AnalysisResult.
    """
    seq = sequence.upper()
    seq_name = name or f"RNA-Pipe-{uuid.uuid4().hex[:6].upper()}"

    # 1.  Composition
    freqs = compute_nucleotide_frequencies(seq)
    composition = NucleotideComposition(
        A=freqs["A"],
        C=freqs["C"],
        G=freqs["G"],
        U=freqs["U"],
    )

    # 2.  Structure prediction
    struct = predict_structure(seq)
    structure = StructurePrediction(
        dot_bracket=struct.dot_bracket,
        minimum_free_energy=struct.minimum_free_energy,
        base_pairs=struct.base_pairs,
        structural_confidence=struct.structural_confidence,
        paired_fraction=struct.paired_fraction,
    )

    # 3.  Stability metrics
    metrics = compute_stability_metrics(
        sequence=seq,
        mfe=struct.minimum_free_energy,
        structural_confidence=struct.structural_confidence,
        paired_fraction=struct.paired_fraction,
        ensemble_diversity=struct.ensemble_diversity,
    )
    stability = StabilityMetrics(
        therapeutic_stability_score=metrics.therapeutic_stability_score,
        gc_content=metrics.gc_content,
        mfe_per_nucleotide=metrics.mfe_per_nucleotide,
        ensemble_diversity=metrics.ensemble_diversity,
    )

    # 4.  Interpretation
    interpretation = interpret_analysis(
        sequence_name=seq_name,
        sequence_length=len(seq),
        gc_content=metrics.gc_content,
        mfe=struct.minimum_free_energy,
        mfe_per_nt=metrics.mfe_per_nucleotide,
        structural_confidence=struct.structural_confidence,
        paired_fraction=struct.paired_fraction,
        ensemble_diversity=metrics.ensemble_diversity,
        therapeutic_stability_score=metrics.therapeutic_stability_score,
        num_base_pairs=len(struct.base_pairs),
    )

    analysis_id = f"ANL-{uuid.uuid4().hex[:8].upper()}"

    return AnalysisResult(
        analysis_id=analysis_id,
        sequence_name=seq_name,
        sequence=seq,
        length=len(seq),
        composition=composition,
        structure=structure,
        stability=stability,
        genomic_interpretation=interpretation,
        timestamp=datetime.utcnow().isoformat() + "Z",
    )


# ---------------------------------------------------------------------------
# 7.  Batch pipeline
# ---------------------------------------------------------------------------


def run_batch_analysis(
    sequences: list[tuple[str, str | None]],
) -> list[AnalysisResult]:
    """Run the full pipeline on a batch of (sequence, name) pairs."""
    return [run_full_analysis(seq, name) for seq, name in sequences]


# ---------------------------------------------------------------------------
# 8.  Export helpers
# ---------------------------------------------------------------------------


def format_fasta(analysis: AnalysisResult) -> str:
    """Format a single analysis as a FASTA entry."""
    desc = (
        f"length={analysis.length} "
        f"GC={analysis.stability.gc_content*100:.1f}% "
        f"stability={analysis.stability.therapeutic_stability_score:.1f}/100 "
        f"MFE={analysis.structure.minimum_free_energy:.2f}"
    )
    return sequence_to_fasta(
        sequence=analysis.sequence,
        name=analysis.sequence_name,
        description=desc,
    )


def format_json_dossier(analyses: list[AnalysisResult]) -> list[dict[str, Any]]:
    """Serialize analyses to JSON-safe dicts."""
    return [a.model_dump() for a in analyses]


# ---------------------------------------------------------------------------
# 9.  Quick aggregate metrics
# ---------------------------------------------------------------------------


def cohort_summary(analyses: list[AnalysisResult]) -> dict[str, Any]:
    """Compute aggregate statistics across a cohort of candidates."""
    scores = [a.stability.therapeutic_stability_score for a in analyses]
    gc_vals = [a.stability.gc_content for a in analyses]
    mfes = [a.structure.minimum_free_energy for a in analyses]
    lengths = [a.length for a in analyses]

    return {
        "count": len(analyses),
        "mean_stability": round(sum(scores) / len(scores), 2) if scores else 0,
        "max_stability": round(max(scores), 2) if scores else 0,
        "min_stability": round(min(scores), 2) if scores else 0,
        "mean_gc": round(sum(gc_vals) / len(gc_vals), 4) if gc_vals else 0,
        "mean_mfe": round(sum(mfes) / len(mfes), 2) if mfes else 0,
        "mean_length": round(sum(lengths) / len(lengths), 1) if lengths else 0,
        "high_viability_count": sum(1 for s in scores if s >= 70),
    }
