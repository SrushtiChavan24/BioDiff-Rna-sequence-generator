"""
RNA Genomics Intelligence Platform — Comparative Analysis Engine.

Generates structured comparative summaries across candidate cohorts,
identifying top performers, strengths, weaknesses, and recommendations
using deterministic multi-metric ranking.
"""

from __future__ import annotations

import math
from typing import Any


def _safe(val: float | None, default: float = 0.0) -> float:
    if val is None:
        return default
    try:
        return val if math.isfinite(val) else default
    except (TypeError, ValueError):
        return default


def generate_comparative_summary(candidates: list[Any]) -> dict:
    """Analyze a cohort of candidates and produce a structured comparison.

    Returns:
        dict with keys:
          - top_candidate: name of the highest-ranked candidate
          - top_analysis_id: analysis ID of the top candidate
          - cohort_size: number of candidates analyzed
          - mean_score: average Therapeutic Stability Index
          - strongest_fold: candidate with best MFE/nt
          - best_gc_balance: candidate with closest-to-optimal GC
          - most_stable_structure: candidate with highest structural confidence
          - highest_therapeutic_score: candidate with max TSI
          - most_flexible_topology: candidate with highest ensemble diversity
          - summary: narrative summary paragraph
          - strengths: list of identified cohort strengths
          - weaknesses: list of identified cohort weaknesses
          - recommendation: professional biotech recommendation
    """
    if not candidates:
        return {
            "top_candidate": "",
            "top_analysis_id": "",
            "cohort_size": 0,
            "mean_score": 0.0,
            "strongest_fold": "",
            "best_gc_balance": "",
            "most_stable_structure": "",
            "highest_therapeutic_score": "",
            "most_flexible_topology": "",
            "summary": "No candidates available for analysis.",
            "strengths": [],
            "weaknesses": [],
            "recommendation": "Submit RNA sequences for analysis to generate comparative insights.",
        }

    scored: list[tuple[float, Any]] = [
        (_safe(c.stability.therapeutic_stability_score), c) for c in candidates
    ]
    scored.sort(key=lambda x: x[0], reverse=True)

    top_candidate = scored[0][1]
    top_score = scored[0][0]
    mean_score = sum(s for s, _ in scored) / len(scored)

    # Strongest fold = best (most negative) MFE per nucleotide
    strongest = min(
        candidates,
        key=lambda c: _safe(c.stability.mfe_per_nucleotide),
    )
    # Best GC balance = closest to 0.50
    best_gc = min(
        candidates,
        key=lambda c: abs(_safe(c.stability.gc_content) - 0.50),
    )
    # Most stable structure = highest structural confidence
    most_confident = max(
        candidates,
        key=lambda c: _safe(c.structure.structural_confidence),
    )
    # Highest therapeutic score
    highest_score = max(
        candidates,
        key=lambda c: _safe(c.stability.therapeutic_stability_score),
    )
    # Most flexible topology = highest ensemble diversity
    most_flexible = max(
        candidates,
        key=lambda c: _safe(c.stability.ensemble_diversity),
    )

    # Identify strengths and weaknesses
    strengths: list[str] = []
    weaknesses: list[str] = []

    avg_gc = sum(_safe(c.stability.gc_content) for c in candidates) / len(candidates)
    avg_conf = sum(_safe(c.structure.structural_confidence) for c in candidates) / len(candidates)
    avg_diversity = sum(_safe(c.stability.ensemble_diversity) for c in candidates) / len(candidates)

    if avg_gc >= 0.45:
        strengths.append(
            f"Cohort GC content ({avg_gc*100:.1f}%) supports thermal stability and nuclease resistance."
        )
    else:
        weaknesses.append(
            f"Cohort GC content ({avg_gc*100:.1f}%) is suboptimal; GC enrichment recommended."
        )

    if mean_score >= 60:
        strengths.append(
            f"Mean Therapeutic Stability Index ({mean_score:.1f}/100) indicates favorable cohort quality."
        )
    elif mean_score >= 40:
        weaknesses.append(
            f"Mean Stability Index ({mean_score:.1f}/100) is moderate; optimization needed."
        )
    else:
        weaknesses.append(
            f"Mean Stability Index ({mean_score:.1f}/100) is low; significant engineering required."
        )

    if avg_conf >= 0.70:
        strengths.append(
            f"High average structural confidence ({avg_conf*100:.1f}%) confirms reliable fold predictions."
        )
    elif avg_conf >= 0.50:
        weaknesses.append(
            f"Moderate structural confidence ({avg_conf*100:.1f}%); some predictions may be ambiguous."
        )
    else:
        weaknesses.append(
            f"Low structural confidence ({avg_conf*100:.1f}%); predictions should be interpreted cautiously."
        )

    if avg_diversity < 0.25:
        strengths.append(
            "Low ensemble diversity across the cohort indicates deterministic folding behavior."
        )
    elif avg_diversity > 0.45:
        weaknesses.append(
            "Elevated ensemble diversity suggests conformational heterogeneity in multiple candidates."
        )

    # Build narrative summary
    top_name = top_candidate.sequence_name
    top_gc = _safe(top_candidate.stability.gc_content) * 100
    top_mfe = _safe(top_candidate.structure.minimum_free_energy)
    top_conf = _safe(top_candidate.structure.structural_confidence) * 100

    summary = (
        f"Comparative analysis of {len(candidates)} RNA candidate(s) identifies "
        f"'{top_name}' as the lead candidate with a Therapeutic Stability Index "
        f"of {top_score:.1f}/100. The lead candidate exhibits a GC content of "
        f"{top_gc:.1f}%, MFE of {top_mfe:.2f} kcal/mol, and structural confidence "
        f"of {top_conf:.1f}%. "
    )

    if strongest.analysis_id != top_candidate.analysis_id:
        summary += (
            f"The strongest thermodynamic fold is observed in "
            f"'{strongest.sequence_name}' (MFE/nt: "
            f"{_safe(strongest.stability.mfe_per_nucleotide):.3f} kcal/mol/nt). "
        )

    if best_gc.analysis_id != top_candidate.analysis_id:
        summary += (
            f"Optimal GC balance is achieved by '{best_gc.sequence_name}' "
            f"({_safe(best_gc.stability.gc_content)*100:.1f}%). "
        )

    if most_confident.analysis_id != top_candidate.analysis_id:
        summary += (
            f"Highest structural reliability belongs to "
            f"'{most_confident.sequence_name}' ({_safe(most_confident.structure.structural_confidence)*100:.1f}% confidence). "
        )

    summary += (
        f"The cohort mean Therapeutic Stability Index is {mean_score:.1f}/100, "
        f"with a mean GC content of {avg_gc*100:.1f}%."
    )

    # Recommendation
    if mean_score >= 60:
        recommendation = (
            f"The candidate cohort demonstrates promising biophysical characteristics. "
            f"Lead candidate '{top_name}' is recommended for advancement to preclinical "
            f"evaluation. Targeted optimization of lower-ranked candidates via GC "
            f"enrichment and stem stabilization is advised to expand the therapeutic pipeline."
        )
    elif mean_score >= 40:
        recommendation = (
            f"The cohort displays moderate therapeutic potential. Lead candidate "
            f"'{top_name}' should undergo sequence optimization focused on enhancing "
            f"conformational rigidity and GC content. Parallel chemical modification "
            f"screening (2'-O-methyl, phosphorothioate) is recommended to identify "
            f"stability improvements."
        )
    else:
        recommendation = (
            f"The current candidate cohort requires substantial biophysical engineering "
            f"before therapeutic application. An iterative design-build-test cycle "
            f"incorporating GC enrichment, stem stabilization, and chemical modification "
            f"is strongly advised. Consider screening a larger sequence library to identify "
            f"candidates with intrinsically favorable folding properties."
        )

    return {
        "top_candidate": top_name,
        "top_analysis_id": top_candidate.analysis_id,
        "cohort_size": len(candidates),
        "mean_score": round(mean_score, 2),
        "strongest_fold": strongest.sequence_name,
        "best_gc_balance": best_gc.sequence_name,
        "most_stable_structure": most_confident.sequence_name,
        "highest_therapeutic_score": highest_score.sequence_name,
        "most_flexible_topology": most_flexible.sequence_name,
        "summary": summary,
        "strengths": strengths,
        "weaknesses": weaknesses,
        "recommendation": recommendation,
    }
