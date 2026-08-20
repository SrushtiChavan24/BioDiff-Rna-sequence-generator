"""
RNA secondary structure analysis utilities.

Pure-function helpers for dot-bracket parsing, base-pair extraction,
and 2D coordinate generation using a naview-inspired folding algorithm.
These are stateless, framework-agnostic, and independently testable.
"""

from __future__ import annotations

import math


def parse_dot_bracket(dot_bracket: str) -> dict:
    """
    Parse dot-bracket notation into a structured representation.

    Returns a dict with:
      - stems: list of stem regions, each containing paired positions
      - loops: list of loop regions (hairpin, internal, bulge)
      - depth_profile: pairing depth at each position
      - base_pairs: list of [i, j] paired indices
      - unpaired_regions: list of [start, end] unpaired segments
    """
    n = len(dot_bracket)
    depth = 0
    depth_profile: list[int] = []
    stack: list[int] = []
    base_pairs: list[list[int]] = []

    for i, ch in enumerate(dot_bracket):
        if ch == "(":
            stack.append(i)
            depth += 1
            depth_profile.append(depth)
        elif ch == ")":
            if stack:
                j = stack.pop()
                base_pairs.append([j, i])
            depth -= 1
            depth_profile.append(max(0, depth))
        else:
            depth_profile.append(depth)

    base_pairs.sort(key=lambda p: p[0])

    stems: list[dict] = []
    loop_regions: list[dict] = []
    unpaired_regions: list[list[int]] = []

    pair_map: dict[int, int] = {}
    for i, j in base_pairs:
        pair_map[i] = j
        pair_map[j] = i

    i = 0
    while i < n:
        if dot_bracket[i] == "(":
            j = pair_map[i]
            stem_pairs: list[list[int]] = []
            k = i
            while k <= j and dot_bracket[k] == "(" and pair_map.get(k, -1) >= k:
                pk = pair_map[k]
                if pk > j:
                    break
                stem_pairs.append([k, pk])
                k += 1
            stems.append({
                "start": i,
                "end": j,
                "pairs": stem_pairs,
                "length": len(stem_pairs),
            })
            last_left = stem_pairs[-1][0]
            first_right = stem_pairs[0][1]
            loop_start = last_left + 1
            loop_end = first_right - 1
            if loop_start <= loop_end:
                loop_regions.append({
                    "type": "hairpin",
                    "start": loop_start,
                    "end": loop_end,
                    "size": loop_end - loop_start + 1,
                })
            i = j + 1
        elif dot_bracket[i] == ".":
            u_start = i
            while i < n and dot_bracket[i] == ".":
                i += 1
            unpaired_regions.append([u_start, i - 1])
        else:
            i += 1

    return {
        "stems": stems,
        "loops": loop_regions,
        "depth_profile": depth_profile,
        "base_pairs": base_pairs,
        "unpaired_regions": unpaired_regions,
        "max_depth": max(depth_profile) if depth_profile else 0,
    }


def extract_base_pairs(dot_bracket: str) -> list[list[int]]:
    """Extract base pair indices from dot-bracket notation using a stack."""
    stack: list[int] = []
    pairs: list[list[int]] = []
    for i, ch in enumerate(dot_bracket):
        if ch == "(":
            stack.append(i)
        elif ch == ")" and stack:
            j = stack.pop()
            pairs.append([j, i])
    pairs.sort(key=lambda p: p[0])
    return pairs


def generate_2d_coordinates(
    sequence: str,
    dot_bracket: str,
    base_pairs: list[list[int]] | None = None,
    width: int = 700,
    height: int = 500,
) -> list[tuple[float, float]]:
    """
    Generate publication-quality 2D coordinates for RNA secondary structure.

    Uses an enhanced naview-inspired radial folding algorithm:
    - Paired nucleotides on parallel stem tracks with natural spacing
    - Loop nucleotides on smooth circular arcs
    - Unpaired regions as organized linear chains
    - Auto-centering of the full structure
    """
    n = len(sequence)
    coords: list[tuple[float, float]] = [(0.0, 0.0)] * n
    placed = [False] * n

    if base_pairs is None:
        base_pairs = extract_base_pairs(dot_bracket)

    pair_map: dict[int, int] = {}
    for i, j in base_pairs:
        pair_map[i] = j
        pair_map[j] = i

    parsed = parse_dot_bracket(dot_bracket)
    stems = parsed["stems"]

    margin = 40.0
    avail_w = width - 2 * margin
    avail_h = height - 2 * margin

    if not stems:
        for i in range(n):
            spacing = min(30.0, max(18.0, 420.0 / max(n, 1)))
            coords[i] = (width / 2 - (n * spacing) / 2 + i * spacing, height / 2 - 20)
            placed[i] = True
        return coords

    n_stems = len(stems)
    stem_spacing_x = avail_w / max(n_stems + 1, 1)

    stem_width = min(70.0, max(50.0, avail_w * 0.12))

    total_height_needed = 0
    for stem in stems:
        sl = len(stem["pairs"])
        total_height_needed = max(total_height_needed, sl)
    vertical_spacing = min(26.0, max(16.0, (avail_h * 0.55) / max(total_height_needed, 1)))

    center_y = height / 2 - 10

    for si, stem in enumerate(stems):
        pairs = stem["pairs"]
        stem_len = len(pairs)

        stem_x = margin + (si + 0.5) * stem_spacing_x

        stem_top_y = center_y - (stem_len * vertical_spacing) / 2

        for k, (li, ri) in enumerate(pairs):
            y_pos = stem_top_y + k * vertical_spacing
            coords[li] = (stem_x, y_pos)
            coords[ri] = (stem_x + stem_width, y_pos)
            placed[li] = True
            placed[ri] = True

        last_left = pairs[-1][0]
        first_right = pairs[0][1]
        loop_start = last_left + 1
        loop_end = first_right - 1

        if loop_start <= loop_end:
            loop_indices = list(range(loop_start, loop_end + 1))
            num_loop = len(loop_indices)
            lx, ly = coords[last_left]
            rx, ry = coords[first_right]

            arc_radius = max(
                32.0,
                (rx - lx) / 2.0 + 14.0,
                num_loop * 16.0 / math.pi,
            )

            for k_idx, idx in enumerate(loop_indices):
                t = (k_idx + 1) / (num_loop + 1)
                angle = -math.pi * 0.9 * t + math.pi * 0.05
                x = (lx + rx) / 2.0 + arc_radius * math.cos(angle + math.pi / 2)
                y = (ly + ry) / 2.0 - arc_radius * 0.05 + arc_radius * math.sin(angle + math.pi / 2) - arc_radius * 0.5
                coords[idx] = (x, y)
                placed[idx] = True

    unplaced = [i for i in range(n) if not placed[i]]
    if unplaced:
        chain_y = height - margin + 10
        chain_spacing = min(24.0, max(14.0, avail_w / max(len(unplaced), 1)))
        chain_start_x = width / 2 - (len(unplaced) * chain_spacing) / 2
        for k, i in enumerate(unplaced):
            coords[i] = (chain_start_x + k * chain_spacing, chain_y)
            placed[i] = True

    xs = [c[0] for c in coords if c[0] != 0 or c[1] != 0 or placed[coords.index(c)]]
    if not xs:
        return coords
    data_min_x = min(x for x, y in coords if placed[coords.index((x, y))])
    data_max_x = max(x for x, y in coords if placed[coords.index((x, y))])
    data_min_y = min(y for x, y in coords if placed[coords.index((x, y))])
    data_max_y = max(y for x, y in coords if placed[coords.index((x, y))])

    data_cx = (data_min_x + data_max_x) / 2
    data_cy = (data_min_y + data_max_y) / 2

    shift_x = width / 2 - data_cx
    shift_y = height / 2 - data_cy - 10

    shifted = []
    for idx, (x, y) in enumerate(coords):
        if placed[idx]:
            shifted.append((x + shift_x, y + shift_y))
        else:
            shifted.append((x, y))

    return shifted


def compute_mountain_profile(dot_bracket: str) -> list[int]:
    """Compute the classic mountain plot profile from dot-bracket notation."""
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


def compute_stem_metrics(dot_bracket: str, base_pairs: list[list[int]]) -> dict:
    """Compute stem-loop structural metrics."""
    parsed = parse_dot_bracket(dot_bracket)
    n = len(dot_bracket)

    stem_count = len(parsed["stems"])
    loop_count = len(parsed["loops"])
    total_stem_pairs = sum(s["length"] for s in parsed["stems"])
    avg_stem_length = total_stem_pairs / max(stem_count, 1)

    loop_sizes = [l["size"] for l in parsed["loops"]]
    avg_loop_size = sum(loop_sizes) / max(len(loop_sizes), 1)
    max_loop_size = max(loop_sizes) if loop_sizes else 0

    paired_pct = (2 * len(base_pairs) / n * 100) if n > 0 else 0

    return {
        "stem_count": stem_count,
        "loop_count": loop_count,
        "total_stem_pairs": total_stem_pairs,
        "avg_stem_length": round(avg_stem_length, 2),
        "avg_loop_size": round(avg_loop_size, 2),
        "max_loop_size": max_loop_size,
        "paired_percentage": round(paired_pct, 1),
        "max_depth": parsed["max_depth"],
    }
