"""
True RNA Secondary Structure Layout and SVG Generation.

Implements a publication-grade RNA secondary structure renderer that
produces authentic stem-loop topology diagrams with curved base-pair
arcs, gradient nucleotide nodes, depth annotations, and professional
biotech aesthetics. Designed for genomics research dashboards and
therapeutic RNA screening systems.
"""

from __future__ import annotations

import math
from xml.sax.saxutils import escape

from app.utils.rna_structure_utils import (
    parse_dot_bracket,
    generate_2d_coordinates,
    compute_stem_metrics,
)

NUC_COLORS = {
    "A": "#4CAF50",
    "U": "#F44336",
    "G": "#FF9800",
    "C": "#2196F3",
}

NUC_GRADIENTS = {
    "A": {"start": "#81C784", "end": "#2E7D32"},
    "U": {"start": "#EF9A9A", "end": "#C62828"},
    "G": {"start": "#FFCC80", "end": "#E65100"},
    "C": {"start": "#90CAF9", "end": "#1565C0"},
}

PAIRED_GLOW_OPACITY = "0.35"
LOOP_FILL_OPACITY = "0.06"
STEM_HIGHLIGHT_OPACITY = "0.08"


def generate_true_fold_svg(
    sequence: str,
    dot_bracket: str,
    base_pairs: list[list[int]],
    width: int = 800,
    height: int = 600,
    title: str | None = None,
) -> str:
    n = len(sequence)
    coords = generate_2d_coordinates(sequence, dot_bracket, base_pairs, width * 2, height * 2)

    pair_map: dict[int, int] = {}
    for i, j in base_pairs:
        pair_map[i] = j
        pair_map[j] = i

    parsed = parse_dot_bracket(dot_bracket)
    depth_profile = parsed["depth_profile"]
    max_depth = parsed["max_depth"] or 1
    stem_metrics = compute_stem_metrics(dot_bracket, base_pairs)

    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    min_x, max_x = (min(xs), max(xs)) if xs else (0, 1)
    min_y, max_y = (min(ys), max(ys)) if ys else (0, 1)

    if max_x - min_x < 1:
        max_x = min_x + 100
    if max_y - min_y < 1:
        max_y = min_y + 100

    pad = 30
    title_h = 36 if title else 0
    legend_h = 38
    avail_w = width - 2 * pad
    avail_h = height - 2 * pad - title_h - legend_h
    data_w = max_x - min_x
    data_h = max_y - min_y

    target_fill = 0.82
    scale = min(
        avail_w * target_fill / max(data_w, 1),
        avail_h * target_fill / max(data_h, 1),
        4.0,
    )
    offset_x = pad + (avail_w - data_w * scale) / 2.0
    offset_y = pad + title_h + (avail_h - data_h * scale) / 2.0 + 4

    def tx(x: float) -> float:
        return offset_x + (x - min_x) * scale

    def ty(y: float) -> float:
        return offset_y + (y - min_y) * scale

    def tx_center(x: float) -> float:
        return offset_x + (x - min_x) * scale

    parts: list[str] = []
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="0 0 {width} {height}" '
        f'width="100%" height="100%" '
        f'style="background-color:#0A0E17;">'
    )

    parts.append("""<defs>
    <filter id="glow" x="-50%" y="-50%" width="200%" height="200%">
        <feGaussianBlur stdDeviation="2.5" result="blur"/>
        <feMerge>
            <feMergeNode in="blur"/>
            <feMergeNode in="SourceGraphic"/>
        </feMerge>
    </filter>
    <filter id="glowStrong" x="-50%" y="-50%" width="200%" height="200%">
        <feGaussianBlur stdDeviation="4" result="blur"/>
        <feMerge>
            <feMergeNode in="blur"/>
            <feMergeNode in="SourceGraphic"/>
        </feMerge>
    </filter>
    <filter id="shadow">
        <feDropShadow dx="0" dy="1" stdDeviation="1.5" flood-color="#000" flood-opacity="0.6"/>
    </filter>
    <linearGradient id="bgGrad" x1="0" y1="0" x2="0" y2="1">
        <stop offset="0%" stop-color="#080B14" stop-opacity="1"/>
        <stop offset="50%" stop-color="#0D1117" stop-opacity="1"/>
        <stop offset="100%" stop-color="#080B14" stop-opacity="1"/>
    </linearGradient>
    <linearGradient id="stemGlow" x1="0" y1="0" x2="0" y2="1">
        <stop offset="0%" stop-color="#7C4DFF" stop-opacity="0"/>
        <stop offset="50%" stop-color="#7C4DFF" stop-opacity="0.08"/>
        <stop offset="100%" stop-color="#7C4DFF" stop-opacity="0"/>
    </linearGradient>
""")

    for nt in "AUGC":
        parts.append(
            f'<radialGradient id="node_{nt}" cx="40%" cy="35%" r="60%">'
            f'<stop offset="0%" stop-color="{NUC_GRADIENTS[nt]["start"]}" stop-opacity="0.95"/>'
            f'<stop offset="70%" stop-color="{NUC_COLORS[nt]}" stop-opacity="0.85"/>'
            f'<stop offset="100%" stop-color="{NUC_GRADIENTS[nt]["end"]}" stop-opacity="0.95"/>'
            f'</radialGradient>'
        )
        parts.append(
            f'<radialGradient id="pairedGlow_{nt}" cx="50%" cy="50%" r="50%">'
            f'<stop offset="0%" stop-color="{NUC_COLORS[nt]}" stop-opacity="0.5"/>'
            f'<stop offset="60%" stop-color="{NUC_COLORS[nt]}" stop-opacity="0.2"/>'
            f'<stop offset="100%" stop-color="{NUC_COLORS[nt]}" stop-opacity="0"/>'
            f'</radialGradient>'
        )

    parts.append("</defs>")

    parts.append(f'<rect x="0" y="0" width="{width}" height="{height}" fill="url(#bgGrad)"/>')

    for stem in parsed["stems"]:
        s_start = stem["pairs"][0][0]
        s_end = stem["pairs"][-1][1]
        if s_start < len(coords) and s_end < len(coords):
            x1 = tx_center(coords[s_start][0]) - 8
            x2 = tx_center(coords[s_end][0]) + 8
            y1 = min(ty(coords[s_start][0]), ty(coords[s_start][1]))
            y2 = max(ty(coords[s_end][0]), ty(coords[s_end][1]))
            if x2 > x1:
                parts.append(
                    f'<rect x="{x1:.1f}" y="{y1 - 4:.1f}" '
                    f'width="{x2 - x1:.1f}" height="{y2 - y1 + 8:.1f}" '
                    f'rx="4" fill="#7C4DFF" fill-opacity="{STEM_HIGHLIGHT_OPACITY}"/>'
                )

    for loop in parsed["loops"]:
        l_start = loop["start"]
        l_end = loop["end"]
        if l_start < len(coords) and l_end < len(coords):
            loop_xs = [tx(coords[i][0]) for i in range(l_start, min(l_end + 1, len(coords)))]
            loop_ys = [ty(coords[i][1]) for i in range(l_start, min(l_end + 1, len(coords)))]
            if loop_xs:
                lx_min, lx_max = min(loop_xs), max(loop_xs)
                ly_min, ly_max = min(loop_ys), max(loop_ys)
                cx = (lx_min + lx_max) / 2
                cy = (ly_min + ly_max) / 2
                rx = (lx_max - lx_min) / 2 + 16
                ry = (ly_max - ly_min) / 2 + 16
                if rx > 0 and ry > 0:
                    parts.append(
                        f'<ellipse cx="{cx:.1f}" cy="{cy:.1f}" '
                        f'rx="{rx:.1f}" ry="{ry:.1f}" '
                        f'fill="#00BCD4" fill-opacity="{LOOP_FILL_OPACITY}" '
                        f'stroke="#00BCD4" stroke-width="0.5" stroke-opacity="0.08" '
                        f'stroke-dasharray="4,3"/>'
                    )

    title_text = title or "RNA Secondary Structure"
    parts.append(
        f'<text x="{width / 2}" y="24" text-anchor="middle" '
        f'fill="#E8EAED" font-family="Inter, \'Segoe UI\', sans-serif" '
        f'font-size="13" font-weight="600" letter-spacing="0.3">'
        f'{escape(title_text)}</text>'
    )

    for i, j in base_pairs:
        if i < j and i < len(coords) and j < len(coords):
            x1, y1 = tx(coords[i][0]), ty(coords[i][1])
            x2, y2 = tx(coords[j][0]), ty(coords[j][1])
            dx = x2 - x1
            dy = y2 - y1
            dist = math.sqrt(dx * dx + dy * dy)
            if dist > 4:
                mid_x = (x1 + x2) / 2
                mid_y = (y1 + y2) / 2
                perp_x = -dy / dist * 18
                perp_y = dx / dist * 18
                cp_x = mid_x + perp_x
                cp_y = mid_y + perp_y - abs(perp_y) * 0.3

                path = (
                    f'M {x1:.1f} {y1:.1f} '
                    f'Q {cp_x:.1f} {cp_y:.1f} '
                    f'{x2:.1f} {y2:.1f}'
                )
                parts.append(
                    f'<path d="{path}" fill="none" '
                    f'stroke="#7C4DFF" stroke-width="0.8" '
                    f'stroke-dasharray="2.5,2" stroke-opacity="0.3" '
                    f'class="base-pair-arc"/>'
                )

    for i in range(n - 1):
        if i < len(coords) and i + 1 < len(coords):
            x1, y1 = tx(coords[i][0]), ty(coords[i][1])
            x2, y2 = tx(coords[i + 1][0]), ty(coords[i + 1][1])
            parts.append(
                f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
                f'stroke="#30363D" stroke-width="1.2" stroke-opacity="0.35" '
                f'class="backbone"/>'
            )

    for i in range(n):
        if i >= len(coords):
            break
        x, y = tx(coords[i][0]), ty(coords[i][1])
        nt = sequence[i] if i < len(sequence) else "?"
        is_paired = i in pair_map
        r = 10 if is_paired else 7
        grad_id = f"node_{nt}" if nt in NUC_GRADIENTS else "node_A"

        if is_paired:
            outer_r = r + 6
            parts.append(
                f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{outer_r}" '
                f'fill="url(#pairedGlow_{nt})" fill-opacity="{PAIRED_GLOW_OPACITY}" '
                f'class="nucleotide-glow" data-index="{i}"/>'
            )

        parts.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" '
            f'fill="url({grad_id})" fill-opacity="0.92" '
            f'stroke="{"#FFFFFF" if is_paired else "#FFFFFF"}" '
            f'stroke-width="{1.2 if is_paired else 0.6}" '
            f'stroke-opacity="{0.35 if is_paired else 0.12}" '
            f'{"filter=\"url(#glow)\"" if is_paired else ""} '
            f'class="nucleotide" data-index="{i}" data-nt="{nt}" '
            f'data-paired="{str(is_paired).lower()}" '
            f'data-depth="{depth_profile[i] if i < len(depth_profile) else 0}"/>'
        )

        font_size = "8.5" if is_paired else "7"
        parts.append(
            f'<text x="{x:.1f}" y="{y + 3}" '
            f'text-anchor="middle" dominant-baseline="central" '
            f'fill="#FFFFFF" font-family="\'JetBrains Mono\', \'Consolas\', monospace" '
            f'font-size="{font_size}" font-weight="700" '
            f'filter="url(#shadow)" '
            f'class="nucleotide-label" data-index="{i}">{nt}</text>'
        )

    legend_y = height - legend_h + 10
    legend_x = pad + 4

    parts.append(
        f'<line x1="{pad}" y1="{legend_y - 10}" '
        f'x2="{width - pad}" y2="{legend_y - 10}" '
        f'stroke="#1A1F2E" stroke-width="0.8"/>'
    )

    parts.append(
        f'<text x="{legend_x}" y="{legend_y}" '
        f'fill="#8B949E" font-family="Inter, sans-serif" '
        f'font-size="7.5" font-weight="500" letter-spacing="1.2">'
        f'LEGEND</text>'
    )
    legend_x += 44

    parts.append(
        f'<circle cx="{legend_x + 5}" cy="{legend_y - 1}" r="5" '
        f'fill="#7C4DFF" fill-opacity="0.2" stroke="#7C4DFF" '
        f'stroke-width="0.8" stroke-opacity="0.4"/>'
    )
    parts.append(
        f'<text x="{legend_x + 14}" y="{legend_y + 2.5}" '
        f'fill="#8B949E" font-family="Inter, sans-serif" '
        f'font-size="8">Paired</text>'
    )
    legend_x += 58
    parts.append(
        f'<circle cx="{legend_x + 5}" cy="{legend_y - 1}" r="4" '
        f'fill="#C9D1D9" fill-opacity="0.2"/>'
    )
    parts.append(
        f'<text x="{legend_x + 13}" y="{legend_y + 2.5}" '
        f'fill="#8B949E" font-family="Inter, sans-serif" '
        f'font-size="8">Unpaired</text>'
    )
    legend_x += 68

    for nt, color in sorted(NUC_COLORS.items()):
        parts.append(
            f'<circle cx="{legend_x + 5}" cy="{legend_y - 1}" r="4.5" '
            f'fill="{color}" fill-opacity="0.85"/>'
        )
        parts.append(
            f'<text x="{legend_x + 13}" y="{legend_y + 2.5}" '
            f'fill="#C9D1D9" font-family="Inter, sans-serif" '
            f'font-size="8" font-weight="500">{nt}</text>'
        )
        legend_x += 28

    parts.append(
        f'<text x="{legend_x + 16}" y="{legend_y + 2.5}" '
        f'fill="#7C4DFF" font-family="Inter, sans-serif" '
        f'font-size="7.5" opacity="0.7">'
        f'Stem: {stem_metrics["stem_count"]} &nbsp;Loop: {stem_metrics["loop_count"]} &nbsp;Depth: {max_depth}</text>'
    )

    num_pairs = len(base_pairs)
    paired_pct = (2 * num_pairs / n * 100) if n > 0 else 0
    stats = (
        f"{n} nt | {num_pairs} bp | "
        f"Paired: {paired_pct:.1f}% | "
        f"Stems: {stem_metrics['stem_count']} | "
        f"Depth: {max_depth}"
    )
    parts.append(
        f'<text x="{width - pad}" y="{legend_y + 2.5}" text-anchor="end" '
        f'fill="#505866" font-family="\'JetBrains Mono\', \'Consolas\', monospace" '
        f'font-size="7.5" opacity="0.8">{escape(stats)}</text>'
    )

    parts.append(
        "<style>"
        ".nucleotide { transition: all 0.12s ease; cursor: pointer; }"
        ".nucleotide:hover { r: 14; stroke: #FFFFFF; stroke-width: 2.5; stroke-opacity: 0.9; filter: url(#glowStrong); }"
        ".nucleotide-glow { transition: opacity 0.12s ease; pointer-events: none; }"
        ".nucleotide-label { pointer-events: none; }"
        ".base-pair-arc { transition: stroke-opacity 0.2s ease; }"
        ".base-pair-arc:hover { stroke-opacity: 0.7; stroke-width: 1.5; }"
        ".backbone { transition: stroke-opacity 0.2s ease; }"
        "</style>"
    )

    parts.append("</svg>")
    return "\n".join(parts)
