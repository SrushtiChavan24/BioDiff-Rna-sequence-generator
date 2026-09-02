"""
RNA Genomics Intelligence Platform — Genomic Interpretation Engine.

Generates concise, computationally grounded interpretations of RNA
secondary structure and thermodynamic stability for therapeutic
candidate screening dossiers.
"""

from __future__ import annotations


def interpret_analysis(
    sequence_name: str,
    sequence_length: int,
    gc_content: float,
    mfe: float,
    mfe_per_nt: float,
    structural_confidence: float,
    paired_fraction: float,
    ensemble_diversity: float,
    therapeutic_stability_score: float,
    num_base_pairs: int,
) -> str:
    """Generate a multi-paragraph genomic interpretation.

    Returns a technically precise, non-repetitive interpretation
    suitable for biotech research documentation and therapeutic
    candidate screening reports.
    """
    paragraphs: list[str] = [
        _sequence_overview(sequence_name, sequence_length, gc_content),
        _structural_assessment(mfe, mfe_per_nt, num_base_pairs, paired_fraction, sequence_length),
        _stability_analysis(therapeutic_stability_score, structural_confidence, ensemble_diversity, gc_content),
    ]
    return "\n\n".join(paragraphs)


def _sequence_overview(name: str, length: int, gc_content: float) -> str:
    gc_pct = gc_content * 100

    if gc_content >= 0.55:
        gc_note = (
            f"elevated GC content ({gc_pct:.1f}%). Triplet hydrogen bonding "
            f"in G-C pairs elevates Tm and confers enhanced nuclease resistance."
        )
    elif gc_content >= 0.40:
        gc_note = (
            f"balanced GC content ({gc_pct:.1f}%) within the optimal therapeutic "
            f"window, supporting stable fold formation while retaining conformational flexibility."
        )
    elif gc_content >= 0.30:
        gc_note = (
            f"moderate GC content ({gc_pct:.1f}%). AU enrichment may reduce "
            f"thermal stability but permits dynamic structural transitions."
        )
    else:
        gc_note = (
            f"low GC content ({gc_pct:.1f}%) indicating AU-rich composition "
            f"associated with reduced duplex stability and potential mRNA lability."
        )

    size_class = _size_classification(length)
    return (
        f"Sequence '{name}' is a {size_class} RNA of {length} nt with {gc_note}"
    )


def _structural_assessment(
    mfe: float,
    mfe_per_nt: float,
    num_base_pairs: int,
    paired_fraction: float,
    length: int,
) -> str:
    paired_pct = paired_fraction * 100

    if mfe_per_nt <= -0.3:
        energy = (
            f"MFE {mfe:.2f} kcal/mol ({mfe_per_nt:.3f} kcal/mol/nt) indicates "
            f"a thermodynamically favorable fold with strong base-pair stacking."
        )
    elif mfe_per_nt <= -0.15:
        energy = (
            f"MFE {mfe:.2f} kcal/mol ({mfe_per_nt:.3f} kcal/mol/nt) reflects "
            f"moderate folding stability with defined stem-loop architecture."
        )
    elif mfe_per_nt <= -0.05:
        energy = (
            f"MFE {mfe:.2f} kcal/mol ({mfe_per_nt:.3f} kcal/mol/nt) suggests "
            f"weak secondary structure with limited base-pairing."
        )
    else:
        energy = (
            f"MFE {mfe:.2f} kcal/mol ({mfe_per_nt:.3f} kcal/mol/nt) indicates "
            f"minimal stable fold formation; predominantly single-stranded."
        )

    return (
        f"{energy} The predicted fold comprises {num_base_pairs} canonical "
        f"base pairs covering {paired_pct:.1f}% of the sequence, defining "
        f"the structural scaffold."
    )


def _stability_analysis(
    score: float,
    confidence: float,
    diversity: float,
    gc_content: float,
) -> str:
    conf_pct = confidence * 100

    if score >= 75:
        tier = "high-viability (Tier I)"
        score_note = (
            f"Therapeutic Stability Index {score:.1f}/100 — "
            f"robust thermodynamics and structural determinacy."
        )
    elif score >= 50:
        tier = "moderate-viability (Tier II)"
        score_note = (
            f"Therapeutic Stability Index {score:.1f}/100 — "
            f"acceptable stability; sequence optimization may improve persistence."
        )
    elif score >= 25:
        tier = "limited-viability (Tier III)"
        score_note = (
            f"Therapeutic Stability Index {score:.1f}/100 — "
            f"marginal stability; structural modifications recommended."
        )
    else:
        tier = "low-viability"
        score_note = (
            f"Therapeutic Stability Index {score:.1f}/100 — "
            f"insufficient stability for direct therapeutic application."
        )

    if confidence >= 0.8:
        conf_note = (
            f"Structural confidence {conf_pct:.1f}% confirms a single dominant fold "
            f"with minimal ensemble heterogeneity (diversity: {diversity:.3f})."
        )
    elif confidence >= 0.5:
        conf_note = (
            f"Structural confidence {conf_pct:.1f}% with moderate ensemble "
            f"diversity ({diversity:.3f}), indicating competing subpopulations."
        )
    else:
        conf_note = (
            f"Structural confidence {conf_pct:.1f}% reflects significant "
            f"conformational heterogeneity (diversity: {diversity:.3f})."
        )

    return f"{score_note} Classified as {tier}. {conf_note}"


def _size_classification(length: int) -> str:
    if length < 50:
        return "short oligonucleotide"
    if length < 200:
        return "medium-length transcript"
    if length < 1000:
        return "long transcript"
    return "extended RNA"
