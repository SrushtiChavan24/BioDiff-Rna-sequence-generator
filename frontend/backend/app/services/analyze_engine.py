from __future__ import annotations

import uuid
from collections import Counter

from app.schemas.analyze_payloads import AnalyzeRequest, AnalyzeResponse
from app.schemas.rna_sequence import RNASequenceInput
from app.services.rna_analysis import analyze_sequence as cache_analysis
from app.services.structure_prediction import predict_structure
from app.services.stability import compute_stability_metrics
from app.interpretations.genomic_interpreter import interpret_analysis
from app.utils.sequence_utils import compute_gc_content


def run_analysis(payload: AnalyzeRequest) -> AnalyzeResponse:
    sequence = payload.sequence
    name = payload.name or f"Candidate-{uuid.uuid4().hex[:6].upper()}"
    n = len(sequence)

    # Persist to analysis cache so structure/candidates pages can retrieve it
    cached = cache_analysis(
        RNASequenceInput(sequence=sequence, name=name, description="via /api/analyze")
    )
    analysis_id = cached.analysis_id

    # 1. Nucleotide counts
    counts = Counter(sequence)
    gc_count = counts.get("G", 0) + counts.get("C", 0)
    au_count = counts.get("A", 0) + counts.get("U", 0)
    gc_content_pct = round((gc_count / n) * 100, 1) if n > 0 else 0.0
    au_content_pct = round((au_count / n) * 100, 1) if n > 0 else 0.0

    # 2. Structure prediction
    struct_result = predict_structure(sequence)
    dot_bracket = struct_result.dot_bracket
    mfe = struct_result.minimum_free_energy
    base_pairs = struct_result.base_pairs
    paired_fraction = struct_result.paired_fraction
    structural_confidence = struct_result.structural_confidence
    ensemble_diversity = struct_result.ensemble_diversity

    # 3. Stability metrics
    gc_ratio = compute_gc_content(sequence)
    stability_result = compute_stability_metrics(
        sequence=sequence,
        mfe=mfe,
        structural_confidence=structural_confidence,
        paired_fraction=paired_fraction,
        ensemble_diversity=ensemble_diversity,
    )

    therapeutic_score = round(stability_result.therapeutic_stability_score, 1)
    structure_confidence_pct = round(structural_confidence * 100, 1)

    # 4. Interpretation
    interpretation = interpret_analysis(
        sequence_name=name,
        sequence_length=n,
        gc_content=gc_ratio,
        mfe=mfe,
        mfe_per_nt=stability_result.mfe_per_nucleotide,
        structural_confidence=structural_confidence,
        paired_fraction=paired_fraction,
        ensemble_diversity=ensemble_diversity,
        therapeutic_stability_score=therapeutic_score,
        num_base_pairs=len(base_pairs),
    )

    candidate_id = f"RNA-{uuid.uuid4().hex[:6].upper()}"

    return AnalyzeResponse(
        candidate_id=candidate_id,
        analysis_id=analysis_id,
        name=name,
        sequence=sequence,
        length=n,
        gc_content=gc_content_pct,
        au_content=au_content_pct,
        mfe=round(mfe, 1),
        paired_ratio=round(paired_fraction, 4),
        therapeutic_score=therapeutic_score,
        structure_confidence=structure_confidence_pct,
        dot_bracket=dot_bracket,
        base_pairs=base_pairs,
        interpretation=interpretation,
    )
