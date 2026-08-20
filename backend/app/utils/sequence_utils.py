"""
RNA sequence analysis utilities.

Pure-function helpers for nucleotide-level sequence characterization.
These are stateless, framework-agnostic, and independently testable.
"""

from __future__ import annotations

from collections import Counter


def compute_nucleotide_frequencies(sequence: str) -> dict[str, float]:
    """
    Compute fractional nucleotide composition.

    Returns dict with keys A, C, G, U and float values summing to 1.0.
    """
    n = len(sequence)
    if n == 0:
        return {"A": 0.0, "C": 0.0, "G": 0.0, "U": 0.0}

    counts = Counter(sequence)
    return {
        "A": counts.get("A", 0) / n,
        "C": counts.get("C", 0) / n,
        "G": counts.get("G", 0) / n,
        "U": counts.get("U", 0) / n,
    }


def compute_gc_content(sequence: str) -> float:
    """
    Compute the guanine-cytosine content ratio.

    GC content is a key thermodynamic indicator: higher GC → stronger
    base pairing → greater thermal stability.
    """
    if not sequence:
        return 0.0
    gc_count = sequence.count("G") + sequence.count("C")
    return gc_count / len(sequence)


def reverse_complement(sequence: str) -> str:
    """
    Compute the reverse complement of an RNA sequence.

    A ↔ U, G ↔ C
    """
    complement_map = str.maketrans("ACGU", "UGCA")
    return sequence.translate(complement_map)[::-1]


def compute_dinucleotide_frequencies(sequence: str) -> dict[str, float]:
    """
    Compute dinucleotide frequencies for sequence complexity analysis.

    Returns fractional counts of all 16 possible dinucleotides.
    """
    n = len(sequence)
    if n < 2:
        return {}

    total = n - 1
    counts: dict[str, int] = {}
    for i in range(total):
        dinuc = sequence[i : i + 2]
        counts[dinuc] = counts.get(dinuc, 0) + 1

    return {k: v / total for k, v in sorted(counts.items())}


def sequence_to_fasta(
    sequence: str,
    name: str = "unnamed_sequence",
    description: str = "",
    line_width: int = 80,
) -> str:
    """
    Format a sequence as a FASTA string.

    Wraps the sequence at the specified line width (default 80).
    """
    header = f">{name}"
    if description:
        header += f" {description}"

    lines = [header]
    for i in range(0, len(sequence), line_width):
        lines.append(sequence[i : i + line_width])

    return "\n".join(lines) + "\n"
