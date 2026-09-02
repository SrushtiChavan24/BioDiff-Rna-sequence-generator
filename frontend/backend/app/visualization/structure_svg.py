"""
RNA Secondary Structure Publication-Style SVG Renderer.

Generates a mountain-plot style SVG visualization of RNA secondary
structure, the standard bioinformatics representation. Shows base-pairing
depth profile across the sequence with color-coded nucleotides below.

This is the canonical RNA secondary structure visualization used in
computational genomics publications (cf. ViennaRNA mountain plots).
"""

from __future__ import annotations

import math
from xml.sax.saxutils import escape

_NUCLEOTIDE_COLORS = {
    "A": "#4CAF50",
    "U": "#F44336",
    "G": "#FF9800",
    "C": "#2196F3",
}

_BG_COLOR = "#0D1117"
_TEXT_COLOR = "#C9D1D9"
_MOUNTAIN_FILL = "#7C4DFF"
_MOUNTAIN_STROKE = "#9C7CFF"
_STEM_HIGHLIGHT = "#00BCD4"
_LOOP_REGION = "#FF6D00"


def _compute_mountain_profile(dot_bracket: str) -> list[int]:
    stack: list[int] = []
    pairs: list[tuple[int, int]] = []
    for i, ch in enumerate(dot_bracket):
        if ch == "(":
            stack.append(i)
        elif ch == ")" and stack:
            j = stack.pop()
            pairs.append((j, i))

    n = len(dot_bracket)
    height = [0] * n
    for j, k in pairs:
        for i in range(j + 1, k):
            height[i] += 1
    return height


def _find_stem_regions(dot_bracket: str) -> list[tuple[int, int, int]]:
    stack: list[int] = []
    stems: list[tuple[int, int, int]] = []
    for i, ch in enumerate(dot_bracket):
        if ch == "(":
            stack.append(i)
        elif ch == ")" and stack:
            j = stack.pop()
            stem_len = i - j + 1
            stems.append((j, i, stem_len))
    return stems


def generate_structure_svg(
    sequence: str,
    dot_bracket: str,
    base_pairs: list[list[int]],
    width: int = 800,
    height: int = 600,
    title: str | None = None,
) -> str:
    n = len(sequence)

    # Layout zones
    top_margin = 60 if title else 40
    mountain_height = max(180, min(300, n // 2))
    seq_height = 80
    bottom_margin = 60
    left_margin = 80
    right_margin = 60
    plot_width = width - left_margin - right_margin

    mountain_top = top_margin
    mountain_bottom = mountain_top + mountain_height
    seq_top = mountain_bottom + 30
    seq_bottom = seq_top + seq_height

    # Compute mountain profile
    profile = _compute_mountain_profile(dot_bracket)
    max_height = max(profile) if profile else 1

    stems = _find_stem_regions(dot_bracket)

    svg_parts: list[str] = []

    svg_parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="0 0 {width} {height}" '
        f'width="{width}" height="{height}" '
        f'style="background-color:{_BG_COLOR};">'
    )

    svg_parts.append("""<defs>
        <linearGradient id="mountGrad" x1="0" y1="0" x2="0" y2="1">
            <stop offset="0%" stop-color="#7C4DFF" stop-opacity="0.35"/>
            <stop offset="100%" stop-color="#7C4DFF" stop-opacity="0.08"/>
        </linearGradient>
        <filter id="glow">
            <feGaussianBlur stdDeviation="1.5" result="blur"/>
            <feMerge><feMergeNode in="blur"/><feMergeNode in="SourceGraphic"/></feMerge>
        </filter>
    </defs>""")

    if title:
        svg_parts.append(
            f'<text x="{width / 2}" y="28" text-anchor="middle" '
            f'fill="{_TEXT_COLOR}" font-family="Inter, sans-serif" '
            f'font-size="15" font-weight="600">{escape(title)}</text>'
        )

    # Y-axis label
    svg_parts.append(
        f'<text x="14" y="{mountain_top + mountain_height / 2}" '
        f'text-anchor="middle" fill="{_TEXT_COLOR}" font-family="Inter, sans-serif" '
        f'font-size="10" transform="rotate(-90, 14, {mountain_top + mountain_height / 2})" '
        f'opacity="0.6">Pairing Depth</text>'
    )

    # Y-axis and gridlines
    y_ticks = 4
    for t in range(y_ticks + 1):
        y_val = mountain_top + mountain_height - (t / y_ticks) * mountain_height
        label_val = int((t / y_ticks) * max_height)
        svg_parts.append(
            f'<text x="{left_margin - 10}" y="{y_val + 3.5}" '
            f'text-anchor="end" fill="{_TEXT_COLOR}" font-family="JetBrains Mono, monospace" '
            f'font-size="8" opacity="0.5">{label_val}</text>'
        )
        if t > 0:
            svg_parts.append(
                f'<line x1="{left_margin}" y1="{y_val}" '
                f'x2="{left_margin + plot_width}" y2="{y_val}" '
                f'stroke="#21262D" stroke-width="0.5" stroke-dasharray="3,3"/>'
            )

    # Mountain fill + stroke
    if n > 0:
        points = []
        for i in range(n):
            x = left_margin + (i / max(1, n - 1)) * plot_width
            h_ratio = profile[i] / max(max_height, 1)
            y = mountain_top + mountain_height - h_ratio * mountain_height
            points.append(f"{x:.1f},{y:.1f}")

        # Fill
        fill_points = (
            f"{left_margin},{mountain_top + mountain_height} "
            + " ".join(points)
            + f" {left_margin + plot_width},{mountain_top + mountain_height}"
        )
        svg_parts.append(
            f'<polygon points="{fill_points}" fill="url(#mountGrad)" stroke="none"/>'
        )

        # Stroke (upper profile)
        svg_parts.append(
            f'<polyline points="{" ".join(points)}" '
            f'fill="none" stroke="{_MOUNTAIN_STROKE}" '
            f'stroke-width="2" stroke-linejoin="round"/>'
        )

        # Stem highlights (vertical markers at stem regions)
        for sj, sk, _ in stems:
            x1 = left_margin + (sj / max(1, n - 1)) * plot_width
            x2 = left_margin + (sk / max(1, n - 1)) * plot_width
            h1 = min(sj, n - 1)
            h2 = min(sk, n - 1)
            y1_p = mountain_top + mountain_height - (profile[h1] / max(max_height, 1)) * mountain_height
            y2_p = mountain_top + mountain_height
            svg_parts.append(
                f'<rect x="{x1}" y="{y1_p}" width="{max(1, x2 - x1)}" '
                f'height="{y2_p - y1_p}" fill="{_STEM_HIGHLIGHT}" '
                f'fill-opacity="0.12" rx="1"/>'
            )
            # Stem bracket at bottom
            svg_parts.append(
                f'<line x1="{x1}" y1="{mountain_top + mountain_height + 4}" '
                f'x2="{x2}" y2="{mountain_top + mountain_height + 4}" '
                f'stroke="{_STEM_HIGHLIGHT}" stroke-width="1.5" stroke-opacity="0.5"/>'
            )
            svg_parts.append(
                f'<line x1="{x1}" y1="{mountain_top + mountain_height + 2}" '
                f'x2="{x1}" y2="{mountain_top + mountain_height + 6}" '
                f'stroke="{_STEM_HIGHLIGHT}" stroke-width="1.5" stroke-opacity="0.5"/>'
            )
            svg_parts.append(
                f'<line x1="{x2}" y1="{mountain_top + mountain_height + 2}" '
                f'x2="{x2}" y2="{mountain_top + mountain_height + 6}" '
                f'stroke="{_STEM_HIGHLIGHT}" stroke-width="1.5" stroke-opacity="0.5"/>'
            )

    # X-axis line
    svg_parts.append(
        f'<line x1="{left_margin}" y1="{mountain_top + mountain_height}" '
        f'x2="{left_margin + plot_width}" y2="{mountain_top + mountain_height}" '
        f'stroke="#30363D" stroke-width="1"/>'
    )

    # X-axis label
    svg_parts.append(
        f'<text x="{left_margin + plot_width / 2}" y="{mountain_top + mountain_height + 16}" '
        f'text-anchor="middle" fill="{_TEXT_COLOR}" font-family="Inter, sans-serif" '
        f'font-size="10" opacity="0.6">Nucleotide Position</text>'
    )

    # Sequence track
    seq_center_y = seq_top + (seq_bottom - seq_top) / 2

    if n <= 200:
        for i in range(n):
            x = left_margin + (i / max(1, n - 1)) * plot_width
            nt = sequence[i] if i < len(sequence) else "?"
            color = _NUCLEOTIDE_COLORS.get(nt, "#999")

            # Determine if paired
            is_paired = any(
                (pi == i or pj == i) for pi, pj in base_pairs
            )

            r = 8 if is_paired else 5
            svg_parts.append(
                f'<circle cx="{x:.1f}" cy="{seq_center_y}" r="{r}" '
                f'fill="{color}" fill-opacity="0.85"'
                f'{" filter=\"url(#glow)\"" if is_paired else ""}/>'
            )
            svg_parts.append(
                f'<text x="{x:.1f}" y="{seq_center_y + 3.5}" '
                f'text-anchor="middle" fill="white" '
                f'font-family="JetBrains Mono, monospace" '
                f'font-size="8" font-weight="bold">{nt}</text>'
            )

            # Connection line to next
            if i < n - 1:
                x2 = left_margin + ((i + 1) / max(1, n - 1)) * plot_width
                svg_parts.append(
                    f'<line x1="{x:.1f}" y1="{seq_center_y}" '
                    f'x2="{x2:.1f}" y2="{seq_center_y}" '
                    f'stroke="#424242" stroke-width="0.8" stroke-opacity="0.3"/>'
                )
    else:
        # Compact mode for long sequences: just show a colored bar
        bar_height = 20
        for i in range(0, n, max(1, n // 200)):
            x = left_margin + (i / max(1, n - 1)) * plot_width
            nt = sequence[i] if i < len(sequence) else "?"
            color = _NUCLEOTIDE_COLORS.get(nt, "#999")
            seg_width = max(2, plot_width / 200)
            svg_parts.append(
                f'<rect x="{x:.1f}" y="{seq_center_y - bar_height / 2}" '
                f'width="{seg_width:.1f}" height="{bar_height}" '
                f'fill="{color}" fill-opacity="0.6" rx="1"/>'
            )

    # Sequence label
    svg_parts.append(
        f'<text x="14" y="{seq_center_y + 3.5}" '
        f'text-anchor="middle" fill="{_TEXT_COLOR}" font-family="Inter, sans-serif" '
        f'font-size="10" opacity="0.6">Sequence</text>'
    )

    # Legend
    legend_y = seq_bottom + 25
    legend_x = left_margin
    for nt, color in _NUCLEOTIDE_COLORS.items():
        svg_parts.append(
            f'<circle cx="{legend_x}" cy="{legend_y}" r="5" fill="{color}" fill-opacity="0.85"/>'
        )
        svg_parts.append(
            f'<text x="{legend_x + 10}" y="{legend_y + 3.5}" '
            f'fill="{_TEXT_COLOR}" font-family="Inter, sans-serif" '
            f'font-size="10">{nt}</text>'
        )
        legend_x += 40

    # Legend: paired/unpaired
    legend_x += 20
    svg_parts.append(
        f'<circle cx="{legend_x}" cy="{legend_y}" r="5" '
        f'fill="{_TEXT_COLOR}" fill-opacity="0.5" filter="url(#glow)"/>'
    )
    svg_parts.append(
        f'<text x="{legend_x + 10}" y="{legend_y + 3.5}" '
        f'fill="{_TEXT_COLOR}" font-family="Inter, sans-serif" '
        f'font-size="10" opacity="0.6">Paired</text>'
    )
    legend_x += 70
    svg_parts.append(
        f'<circle cx="{legend_x}" cy="{legend_y}" r="4" '
        f'fill="{_TEXT_COLOR}" fill-opacity="0.3"/>'
    )
    svg_parts.append(
        f'<text x="{legend_x + 10}" y="{legend_y + 3.5}" '
        f'fill="{_TEXT_COLOR}" font-family="Inter, sans-serif" '
        f'font-size="10" opacity="0.6">Unpaired</text>'
    )

    # Stats footer
    num_pairs = len(base_pairs)
    paired_pct = (2 * num_pairs / n * 100) if n > 0 else 0
    stats = (
        f"{n} nt | {num_pairs} bp | "
        f"Paired: {paired_pct:.1f}% | "
        f"Mountain max: {max_height}"
    )
    svg_parts.append(
        f'<text x="{width - right_margin}" y="{height - 15}" text-anchor="end" '
        f'fill="{_TEXT_COLOR}" font-family="JetBrains Mono, monospace" '
        f'font-size="9" opacity="0.5">{escape(stats)}</text>'
    )

    svg_parts.append("</svg>")
    return "\n".join(svg_parts)
