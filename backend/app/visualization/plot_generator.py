"""
Plotly Chart Specification Generator.

Generates Plotly JSON chart definitions server-side for frontend
rendering. All charts use a consistent dark genomics theme.
"""

from __future__ import annotations

import numpy as np


# Consistent dark theme for all charts
_DARK_THEME = {
    "paper_bgcolor": "#0D1117",
    "plot_bgcolor": "#161B22",
    "font": {"color": "#C9D1D9", "family": "Inter, sans-serif"},
    "xaxis": {
        "gridcolor": "#21262D",
        "linecolor": "#30363D",
        "zerolinecolor": "#30363D",
    },
    "yaxis": {
        "gridcolor": "#21262D",
        "linecolor": "#30363D",
        "zerolinecolor": "#30363D",
    },
}

# Biotech color palette
_FAMILY_COLORS = [
    "#7C4DFF", "#00BCD4", "#FF6D00", "#64DD17",
    "#FF4081", "#FFD600", "#00E5FF", "#E040FB",
    "#76FF03", "#FF3D00",
]


def generate_similarity_heatmap(
    matrix: np.ndarray,
    labels: list[str],
) -> dict:
    """
    Generate Plotly spec for a molecular similarity landscape heatmap.

    Returns a Plotly figure dictionary ready for JSON serialization.
    """
    return {
        "data": [
            {
                "type": "heatmap",
                "z": matrix.tolist(),
                "x": labels,
                "y": labels,
                "colorscale": [
                    [0.0, "#0D1117"],
                    [0.25, "#1A237E"],
                    [0.5, "#7C4DFF"],
                    [0.75, "#B388FF"],
                    [1.0, "#E8EAF6"],
                ],
                "colorbar": {
                    "title": {"text": "Similarity", "font": {"size": 12}},
                    "tickfont": {"size": 10},
                },
                "hovertemplate": (
                    "<b>%{x}</b> vs <b>%{y}</b><br>"
                    "Similarity: %{z:.3f}<extra></extra>"
                ),
            }
        ],
        "layout": {
            **_DARK_THEME,
            "title": {
                "text": "Molecular Similarity Landscape",
                "font": {"size": 18},
            },
            "width": 700,
            "height": 600,
            "xaxis": {
                **_DARK_THEME["xaxis"],
                "tickangle": -45,
                "tickfont": {"size": 10},
            },
            "yaxis": {
                **_DARK_THEME["yaxis"],
                "tickfont": {"size": 10},
                "autorange": "reversed",
            },
            "margin": {"l": 120, "r": 40, "t": 60, "b": 120},
        },
    }


def generate_projection_scatter(
    points: list[dict],
    family_labels: dict[int, str],
) -> dict:
    """
    Generate Plotly spec for a structural relationship projection scatter.

    Args:
        points: List of dicts with x, y, family_id, label keys.
        family_labels: Map from family_id to family name.
    """
    # Group by family
    families: dict[int, list[dict]] = {}
    for pt in points:
        fid = pt["family_id"]
        families.setdefault(fid, []).append(pt)

    traces = []
    for fid, members in sorted(families.items()):
        color = _FAMILY_COLORS[fid % len(_FAMILY_COLORS)]
        label = family_labels.get(fid, f"Family {fid}")

        traces.append(
            {
                "type": "scatter",
                "mode": "markers+text",
                "x": [p["x"] for p in members],
                "y": [p["y"] for p in members],
                "text": [p["label"] for p in members],
                "textposition": "top center",
                "textfont": {"size": 9, "color": color},
                "name": label,
                "marker": {
                    "size": 12,
                    "color": color,
                    "opacity": 0.85,
                    "line": {"width": 1, "color": "#ffffff30"},
                },
                "hovertemplate": (
                    f"<b>%{{text}}</b><br>"
                    f"{label}<br>"
                    f"X: %{{x:.3f}}<br>"
                    f"Y: %{{y:.3f}}<extra></extra>"
                ),
            }
        )

    return {
        "data": traces,
        "layout": {
            **_DARK_THEME,
            "title": {
                "text": "Structural Relationship Projection",
                "font": {"size": 18},
            },
            "xaxis": {
                **_DARK_THEME["xaxis"],
                "title": {"text": "Component 1"},
            },
            "yaxis": {
                **_DARK_THEME["yaxis"],
                "title": {"text": "Component 2"},
            },
            "width": 700,
            "height": 550,
            "showlegend": True,
            "legend": {
                "font": {"size": 11},
                "bgcolor": "#161B22",
                "bordercolor": "#30363D",
                "borderwidth": 1,
            },
            "margin": {"l": 60, "r": 40, "t": 60, "b": 60},
        },
    }


def generate_stability_comparison(
    names: list[str],
    scores: list[float],
    gc_contents: list[float],
    mfe_values: list[float],
) -> dict:
    """Generate a grouped bar chart comparing stability metrics across candidates."""
    return {
        "data": [
            {
                "type": "bar",
                "name": "Therapeutic Stability Score",
                "x": names,
                "y": scores,
                "marker": {"color": "#7C4DFF", "opacity": 0.85},
                "hovertemplate": "<b>%{x}</b><br>Score: %{y:.1f}/100<extra></extra>",
            },
            {
                "type": "bar",
                "name": "GC Content (%)",
                "x": names,
                "y": [gc * 100 for gc in gc_contents],
                "marker": {"color": "#00BCD4", "opacity": 0.85},
                "hovertemplate": "<b>%{x}</b><br>GC: %{y:.1f}%<extra></extra>",
            },
        ],
        "layout": {
            **_DARK_THEME,
            "title": {
                "text": "Therapeutic Candidate Comparison",
                "font": {"size": 18},
            },
            "barmode": "group",
            "xaxis": {
                **_DARK_THEME["xaxis"],
                "tickangle": -30,
            },
            "yaxis": {
                **_DARK_THEME["yaxis"],
                "title": {"text": "Value"},
            },
            "width": 700,
            "height": 450,
            "showlegend": True,
            "legend": {
                "font": {"size": 11},
                "bgcolor": "#161B22",
            },
            "margin": {"l": 60, "r": 40, "t": 60, "b": 80},
        },
    }


def generate_nucleotide_composition_chart(
    names: list[str],
    compositions: list[dict[str, float]],
) -> dict:
    """Generate stacked bar chart of nucleotide composition across sequences."""
    nucleotides = ["A", "C", "G", "U"]
    colors = {"A": "#4CAF50", "C": "#2196F3", "G": "#FF9800", "U": "#F44336"}

    traces = []
    for nt in nucleotides:
        traces.append(
            {
                "type": "bar",
                "name": nt,
                "x": names,
                "y": [comp.get(nt, 0) * 100 for comp in compositions],
                "marker": {"color": colors[nt], "opacity": 0.85},
                "hovertemplate": f"<b>%{{x}}</b><br>{nt}: %{{y:.1f}}%<extra></extra>",
            }
        )

    return {
        "data": traces,
        "layout": {
            **_DARK_THEME,
            "title": {
                "text": "Genomic Feature Representation",
                "font": {"size": 18},
            },
            "barmode": "stack",
            "xaxis": {
                **_DARK_THEME["xaxis"],
                "tickangle": -30,
            },
            "yaxis": {
                **_DARK_THEME["yaxis"],
                "title": {"text": "Composition (%)"},
                "range": [0, 100],
            },
            "width": 700,
            "height": 450,
            "showlegend": True,
            "legend": {"font": {"size": 11}, "bgcolor": "#161B22"},
            "margin": {"l": 60, "r": 40, "t": 60, "b": 80},
        },
    }
