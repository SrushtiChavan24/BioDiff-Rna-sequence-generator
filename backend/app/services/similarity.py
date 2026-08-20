"""
Molecular Similarity Computation Service.

Computes pairwise similarity between RNA sequences using a multi-metric
approach combining:
- Sequence-level edit distance (Levenshtein normalized)
- Structural distance (Hamming on dot-bracket notation)
- Compositional similarity (nucleotide frequency correlation)

The result is an N×N similarity matrix where each entry is in [0, 1]
with 1.0 indicating identical sequences.
"""

from __future__ import annotations

import numpy as np

from app.utils.sequence_utils import compute_nucleotide_frequencies


def _normalized_edit_distance(seq1: str, seq2: str) -> float:
    """
    Compute normalized Levenshtein edit distance between two sequences.

    Returns a value in [0, 1] where 0 = identical, 1 = completely different.
    Uses O(min(m,n)) space optimization.
    """
    m, n = len(seq1), len(seq2)
    if m == 0 and n == 0:
        return 0.0
    if m == 0 or n == 0:
        return 1.0

    # Ensure seq1 is the shorter one for space optimization
    if m > n:
        seq1, seq2 = seq2, seq1
        m, n = n, m

    prev = list(range(m + 1))
    curr = [0] * (m + 1)

    for j in range(1, n + 1):
        curr[0] = j
        for i in range(1, m + 1):
            cost = 0 if seq1[i - 1] == seq2[j - 1] else 1
            curr[i] = min(
                curr[i - 1] + 1,      # insertion
                prev[i] + 1,          # deletion
                prev[i - 1] + cost,   # substitution
            )
        prev, curr = curr, prev

    return prev[m] / max(m, n)


def _structural_distance(struct1: str, struct2: str) -> float:
    """
    Compute structural distance between two dot-bracket strings.

    Uses padded Hamming distance normalized by the longer structure.
    """
    m, n = len(struct1), len(struct2)
    max_len = max(m, n)
    if max_len == 0:
        return 0.0

    # Pad shorter structure with dots (unpaired)
    s1 = struct1.ljust(max_len, ".")
    s2 = struct2.ljust(max_len, ".")

    mismatches = sum(1 for a, b in zip(s1, s2) if a != b)
    return mismatches / max_len


def _compositional_similarity(seq1: str, seq2: str) -> float:
    """
    Compute compositional similarity based on nucleotide frequencies.

    Uses cosine similarity between frequency vectors.
    """
    freq1 = compute_nucleotide_frequencies(seq1)
    freq2 = compute_nucleotide_frequencies(seq2)

    vec1 = np.array([freq1["A"], freq1["C"], freq1["G"], freq1["U"]])
    vec2 = np.array([freq2["A"], freq2["C"], freq2["G"], freq2["U"]])

    dot = np.dot(vec1, vec2)
    norm1 = np.linalg.norm(vec1)
    norm2 = np.linalg.norm(vec2)

    if norm1 == 0 or norm2 == 0:
        return 0.0

    return float(dot / (norm1 * norm2))


def compute_similarity_matrix(
    sequences: list[str],
    structures: list[str],
    weights: dict[str, float] | None = None,
) -> np.ndarray:
    """
    Compute an N×N molecular similarity matrix.

    Args:
        sequences: List of RNA sequences.
        structures: Corresponding dot-bracket structures.
        weights: Optional dict with keys 'sequence', 'structure', 'composition'.
                 Defaults to equal weights.

    Returns:
        N×N numpy array with similarity values in [0, 1].
    """
    if weights is None:
        weights = {"sequence": 0.4, "structure": 0.4, "composition": 0.2}

    n = len(sequences)
    matrix = np.eye(n, dtype=np.float64)  # Diagonal = 1.0 (self-similarity)

    for i in range(n):
        for j in range(i + 1, n):
            # Sequence similarity (1 - edit distance)
            seq_sim = 1.0 - _normalized_edit_distance(sequences[i], sequences[j])

            # Structural similarity (1 - structural distance)
            struct_sim = 1.0 - _structural_distance(structures[i], structures[j])

            # Compositional similarity
            comp_sim = _compositional_similarity(sequences[i], sequences[j])

            # Weighted composite
            similarity = (
                weights["sequence"] * seq_sim
                + weights["structure"] * struct_sim
                + weights["composition"] * comp_sim
            )

            matrix[i][j] = round(similarity, 4)
            matrix[j][i] = round(similarity, 4)

    return matrix
