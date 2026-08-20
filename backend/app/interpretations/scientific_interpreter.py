"""
RNA Genomics Intelligence Platform — Advanced Scientific Interpretation Engine.

Generates deterministic, rule-based biotech-grade scientific interpretations
from computed RNA metrics. No AI/LLM dependency — all logic is algorithmic.
"""

from __future__ import annotations

import math
from typing import Any

from app.utils.rna_structure_utils import compute_stem_metrics


def _safe(val: float | None, default: float = 0.0) -> float:
    if val is None:
        return default
    try:
        return val if math.isfinite(val) else default
    except (TypeError, ValueError):
        return default


def generate_scientific_interpretation(candidate: Any) -> str:
    """Generate a comprehensive biotech-grade scientific interpretation.

    Returns a multi-paragraph deterministic interpretation covering GC balance,
    structural stability, thermodynamic viability, translational persistence,
    folding compactness, ensemble diversity, stem-loop complexity,
    therapeutic suitability, structural confidence, and degradation resistance.
    """
    seq = candidate.sequence or ""
    seq_len = _safe(float(candidate.length)) if hasattr(candidate, "length") else len(seq)

    gc = _safe(candidate.stability.gc_content)
    gc_pct = gc * 100
    score = _safe(candidate.stability.therapeutic_stability_score)
    mfe = _safe(candidate.structure.minimum_free_energy)
    mfe_per_nt = _safe(candidate.stability.mfe_per_nucleotide)
    conf = _safe(candidate.structure.structural_confidence)
    conf_pct = conf * 100
    paired = _safe(candidate.structure.paired_fraction)
    paired_pct = paired * 100
    diversity = _safe(candidate.stability.ensemble_diversity)
    num_pairs = len(candidate.structure.base_pairs or [])

    dot_bracket = candidate.structure.dot_bracket or ""
    base_pairs = candidate.structure.base_pairs or []
    stem_metrics = compute_stem_metrics(dot_bracket, base_pairs)

    paragraphs: list[str] = []

    # 1. GC balance assessment
    if gc >= 0.55:
        gc_assessment = (
            f"GC Content Assessment: The sequence exhibits elevated guanine-cytosine "
            f"enrichment ({gc_pct:.1f}%), characteristic of thermodynamically resilient "
            f"RNA structures. The triplet hydrogen bonding in G-C pairs elevates the "
            f"melting temperature (Tm) and confers enhanced resistance to endonuclease "
            f"degradation, a desirable property for therapeutic RNA candidates requiring "
            f"prolonged intracellular persistence."
        )
    elif gc >= 0.40:
        gc_assessment = (
            f"GC Content Assessment: The sequence displays a balanced GC composition "
            f"({gc_pct:.1f}%), falling within the optimal window for therapeutic RNA "
            f"applications. This ratio supports stable secondary structure formation "
            f"while preserving the conformational plasticity necessary for target "
            f"recognition and protein interaction."
        )
    elif gc >= 0.30:
        gc_assessment = (
            f"GC Content Assessment: The GC fraction ({gc_pct:.1f}%) is moderately "
            f"below the optimal therapeutic range. The corresponding AU enrichment "
            f"facilitates dynamic structural transitions but may reduce baseline "
            f"thermal stability and increase susceptibility to RNase-mediated decay."
        )
    else:
        gc_assessment = (
            f"GC Content Assessment: The sequence is AU-rich with a low GC fraction "
            f"({gc_pct:.1f}%). This composition predisposes the transcript to reduced "
            f"duplex stability and accelerated turnover via AU-rich element (ARE) "
            f"mediated decay pathways. Chemical modification strategies (2'-O-methyl, "
            f"phosphorothioate linkages) are strongly advised for therapeutic applications."
        )
    paragraphs.append(gc_assessment)

    # 2. Structural stability & thermodynamic viability
    if mfe_per_nt <= -0.30:
        thermo_assessment = (
            f"Thermodynamic Stability: The minimum free energy of {mfe:.2f} kcal/mol "
            f"(normalized: {mfe_per_nt:.3f} kcal/mol/nt) signifies a strongly favored "
            f"folding landscape with deep thermodynamic minima. The base-pairing density "
            f"({paired_pct:.1f}%) confirms extensive stacking interactions and a compact, "
            f"well-ordered tertiary fold. This energetic profile is consistent with "
            f"conformational rigidity under physiological ionic conditions."
        )
    elif mfe_per_nt <= -0.15:
        thermo_assessment = (
            f"Thermodynamic Stability: The predicted fold carries a minimum free energy "
            f"of {mfe:.2f} kcal/mol ({mfe_per_nt:.3f} kcal/mol/nt), indicative of "
            f"moderate thermodynamic stabilization. The structure achieves "
            f"{paired_pct:.1f}% pairing coverage across {num_pairs} base pairs, "
            f"forming a recognizable stem-loop architecture. Competing low-energy "
            f"conformations may coexist within 1-2 kcal/mol of the ground state."
        )
    elif mfe_per_nt <= -0.05:
        thermo_assessment = (
            f"Thermodynamic Stability: The MFE of {mfe:.2f} kcal/mol "
            f"({mfe_per_nt:.3f} kcal/mol/nt) suggests a weakly structured "
            f"conformation with limited base-pairing ({paired_pct:.1f}% paired). "
            f"The shallow energy landscape permits conformational sampling across "
            f"multiple near-ground-state structures, which may influence target "
            f"accessibility and biological activity."
        )
    else:
        thermo_assessment = (
            f"Thermodynamic Stability: The MFE of {mfe:.2f} kcal/mol "
            f"({mfe_per_nt:.3f} kcal/mol/nt) indicates minimal thermodynamic "
            f"driving force for secondary structure formation. The predominant "
            f"single-stranded conformation maximizes sequence availability for "
            f"hybridization but offers limited structural protection against "
            f"enzymatic degradation."
        )
    paragraphs.append(thermo_assessment)

    # 3. Translational persistence & folding compactness
    if paired_pct >= 60:
        compactness = (
            f"Structural Compactness: With {paired_pct:.1f}% of nucleotides "
            f"participating in base pairing, the RNA adopts a densely packed "
            f"topology. The extensive hydrogen-bond network creates a structurally "
            f"rigid scaffold that may impede ribosomal scanning if located within "
            f"the 5' untranslated region. The competing demands of structural "
            f"stability and translational accessibility should be considered in "
            f"therapeutic design."
        )
    elif paired_pct >= 35:
        compactness = (
            f"Structural Compactness: A pairing fraction of {paired_pct:.1f}% "
            f"demarcates {stem_metrics['stem_count']} stem region(s) interspersed "
            f"with {stem_metrics['loop_count']} loop structure(s). This balanced "
            f"architecture provides structural scaffolding for intracellular "
            f"persistence while preserving single-stranded segments for "
            f"ribosomal engagement and RNA-protein interactions."
        )
    else:
        compactness = (
            f"Structural Compactness: The limited pairing density ({paired_pct:.1f}%) "
            f"reflects an open, minimally structured conformation conducive to "
            f"efficient translation initiation. However, the absence of extensive "
            f"secondary structure reduces the RNA's resistance to exonuclease "
            f"activity and may compromise intracellular half-life."
        )
    paragraphs.append(compactness)

    # 4. Ensemble diversity & structural determinacy
    if diversity < 0.15:
        diversity_assessment = (
            f"Ensemble Determinacy: The ensemble diversity metric ({diversity:.3f}) "
            f"confirms a highly deterministic folding pathway with a single "
            f"dominant conformation comprising >90% of the Boltzmann-weighted "
            f"ensemble. This fold homogeneity is strongly favorable for therapeutic "
            f"predictability, ensuring consistent target engagement and uniform "
            f"pharmacokinetic behavior across dosing iterations."
        )
    elif diversity < 0.35:
        diversity_assessment = (
            f"Ensemble Determinacy: Moderate ensemble diversity ({diversity:.3f}) "
            f"indicates a primary fold population with measurable subpopulations "
            f"of alternative conformations. While the dominant structure accounts "
            f"for the majority of the ensemble, conformational switching between "
            f"near-isoenergetic states may introduce variability in target "
            f"accessibility."
        )
    else:
        diversity_assessment = (
            f"Ensemble Determinacy: The elevated ensemble diversity ({diversity:.3f}) "
            f"reveals significant conformational heterogeneity with multiple "
            f"competing structural states. This thermodynamic polymorphism may "
            f"lead to unpredictable in vivo folding behavior, potentially "
            f"compromising both efficacy and safety profiles in therapeutic contexts."
        )
    paragraphs.append(diversity_assessment)

    # 5. Therapeutic suitability & degradation resistance
    if score >= 70:
        therapeutic = (
            f"Therapeutic Suitability: The composite Therapeutic Stability Index "
            f"({score:.1f}/100) classifies this candidate as a high-priority lead. "
            f"Structural confidence ({conf_pct:.1f}%) supports reliable fold "
            f"prediction, and the stem-loop topology "
            f"({stem_metrics['stem_count']} stems, {stem_metrics['loop_count']} loops) "
            f"provides a robust structural platform. The candidate is recommended "
            f"for advancement to preclinical validation studies, with attention to "
            f"delivery formulation optimization."
        )
    elif score >= 45:
        therapeutic = (
            f"Therapeutic Suitability: With a Therapeutic Stability Index of "
            f"{score:.1f}/100, this candidate falls within the moderate viability "
            f"range. Structural confidence ({conf_pct:.1f}%) provides reasonable "
            f"assurance in the predicted fold ({stem_metrics['stem_count']} stems, "
            f"{stem_metrics['loop_count']} loops). Targeted sequence optimization "
            f"directed at regions of low ensemble stability is recommended to "
            f"enhance conformational rigidity and nuclease resistance."
        )
    else:
        therapeutic = (
            f"Therapeutic Suitability: The Therapeutic Stability Index "
            f"({score:.1f}/100) indicates suboptimal biophysical characteristics "
            f"for direct therapeutic application. Structural confidence "
            f"({conf_pct:.1f}%) suggests the predicted fold should be interpreted "
            f"with caution. Iterative sequence redesign incorporating GC "
            f"enrichment and stem stabilization motifs, combined with chemical "
            f"modification (2'-O-methyl, phosphorothioate), is advised to achieve "
            f"therapeutic-grade stability."
        )
    paragraphs.append(therapeutic)

    # 6. Degradation resistance estimation
    if gc_pct >= 55 and paired_pct >= 50:
        degradation = (
            f"Degradation Resistance: The combination of elevated GC content "
            f"({gc_pct:.1f}%) and extensive base pairing ({paired_pct:.1f}%) "
            f"confers strong resistance to both endonuclease and exonuclease "
            f"activity. The predicted intracellular half-life is expected to "
            f"exceed that of typical unmodified RNA transcripts, supporting "
            f"sustained pharmacological activity."
        )
    elif gc_pct >= 40 and paired_pct >= 35:
        degradation = (
            f"Degradation Resistance: Moderate protection against nuclease-mediated "
            f"decay is expected given the GC content ({gc_pct:.1f}%) and pairing "
            f"density ({paired_pct:.1f}%). Incorporation of modified nucleotides "
            f"at exposed loop regions is recommended to achieve the degradation "
            f"resistance profile required for clinical development."
        )
    else:
        degradation = (
            f"Degradation Resistance: The sequence exhibits limited intrinsic "
            f"resistance to nuclease degradation, attributable to the lower GC "
            f"content ({gc_pct:.1f}%) and sparse secondary structure "
            f"({paired_pct:.1f}% paired). A comprehensive chemical modification "
            f"strategy encompassing backbone and sugar modifications is essential "
            f"for therapeutic application."
        )
    paragraphs.append(degradation)

    return "\n\n".join(paragraphs)
