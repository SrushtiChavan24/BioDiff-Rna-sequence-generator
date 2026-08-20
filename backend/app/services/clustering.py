"""
Structural Family Detection Service.

Identifies candidate structural families among analyzed RNA sequences
using agglomerative hierarchical grouping on the molecular similarity matrix,
and projects sequences into a 2D structural relationship space.

Terminology mapping:
- "Clusters" → "Candidate Structural Families"
- "PCA" → "Structural Relationship Projection"
- "Latent Space" → "Sequence Relationship Space"
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


@dataclass
class FamilyInfo:
    """Internal representation of a structural family."""

    family_id: int
    family_label: str
    member_indices: list[int]
    centroid_index: int
    avg_similarity: float


@dataclass
class ProjectionPointData:
    """A point in the 2D structural relationship projection."""

    index: int
    x: float
    y: float
    family_id: int


@dataclass
class ClusteringOutput:
    """Complete output of structural family detection."""

    families: list[FamilyInfo]
    projection: list[ProjectionPointData]
    total_families: int


def _determine_num_families(n_samples: int) -> int:
    """
    Heuristically determine the number of structural families.

    For small datasets: fewer families.
    For larger: scale logarithmically up to a cap.
    """
    if n_samples <= 2:
        return 1
    if n_samples <= 5:
        return min(2, n_samples)
    if n_samples <= 10:
        return min(3, n_samples)
    return min(int(np.log2(n_samples)) + 1, 10)


def detect_structural_families(
    similarity_matrix: np.ndarray,
    labels: list[str],
    n_families: int | None = None,
) -> ClusteringOutput:
    """
    Detect candidate structural families from a similarity matrix.

    Args:
        similarity_matrix: N×N similarity matrix (values 0–1).
        labels: Sequence labels for identification.
        n_families: Number of families to detect. Auto-determined if None.

    Returns:
        ClusteringOutput with family assignments and 2D projection.
    """
    n = len(labels)

    if n <= 1:
        return ClusteringOutput(
            families=[
                FamilyInfo(
                    family_id=0,
                    family_label="Structural Family Alpha",
                    member_indices=[0] if n == 1 else [],
                    centroid_index=0,
                    avg_similarity=1.0,
                )
            ],
            projection=[
                ProjectionPointData(index=0, x=0.0, y=0.0, family_id=0)
            ]
            if n == 1
            else [],
            total_families=1,
        )

    # Convert similarity to distance for clustering
    distance_matrix = 1.0 - similarity_matrix
    np.fill_diagonal(distance_matrix, 0.0)

    # Determine number of families
    if n_families is None:
        n_families = _determine_num_families(n)
    n_families = min(n_families, n)

    # Agglomerative clustering with precomputed distances
    clustering = AgglomerativeClustering(
        n_clusters=n_families,
        metric="precomputed",
        linkage="average",
    )
    cluster_labels = clustering.fit_predict(distance_matrix)

    # Greek-letter family naming
    family_names = [
        "Alpha", "Beta", "Gamma", "Delta", "Epsilon",
        "Zeta", "Eta", "Theta", "Iota", "Kappa",
    ]

    # Build family info
    families: list[FamilyInfo] = []
    for fid in range(n_families):
        members = [i for i, cl in enumerate(cluster_labels) if cl == fid]
        if not members:
            continue

        # Find centroid (member with highest average similarity to others)
        if len(members) == 1:
            centroid = members[0]
            avg_sim = 1.0
        else:
            sub_matrix = similarity_matrix[np.ix_(members, members)]
            avg_sims = sub_matrix.mean(axis=1)
            centroid_local = int(np.argmax(avg_sims))
            centroid = members[centroid_local]
            # Average pairwise similarity within family
            upper_tri = sub_matrix[np.triu_indices(len(members), k=1)]
            avg_sim = float(upper_tri.mean()) if len(upper_tri) > 0 else 1.0

        name = family_names[fid] if fid < len(family_names) else f"Family-{fid + 1}"

        families.append(
            FamilyInfo(
                family_id=fid,
                family_label=f"Structural Family {name}",
                member_indices=members,
                centroid_index=centroid,
                avg_similarity=round(avg_sim, 4),
            )
        )

    # 2D Structural Relationship Projection (therapeutic topology mapping)
    scaler = StandardScaler()
    scaled = scaler.fit_transform(similarity_matrix)

    n_components = min(2, n, similarity_matrix.shape[1])
    pca = PCA(n_components=n_components)
    coords = pca.fit_transform(scaled)

    projection: list[ProjectionPointData] = []
    for i in range(n):
        projection.append(
            ProjectionPointData(
                index=i,
                x=round(float(coords[i, 0]), 4),
                y=round(float(coords[i, 1]) if n_components > 1 else 0.0, 4),
                family_id=int(cluster_labels[i]),
            )
        )

    return ClusteringOutput(
        families=families,
        projection=projection,
        total_families=len(families),
    )
