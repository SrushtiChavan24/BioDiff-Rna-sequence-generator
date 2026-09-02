"""
Biotech-Grade RNA Therapeutic Dossier PDF Generator.

Produces a professional genomic intelligence report with dark scientific
styling, embedded structure visualizations, analytical charts, and
comprehensive therapeutic candidate assessment.
"""

from __future__ import annotations

import io
import math
from datetime import datetime
from typing import Any

from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
    PageBreak,
    Flowable,
    HRFlowable,
)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.pdfgen import canvas

from app.schemas.analysis_result import AnalysisResult
from app.utils.rna_structure_utils import (
    generate_2d_coordinates,
    compute_mountain_profile,
    compute_stem_metrics,
)
from app.interpretations.scientific_interpreter import generate_scientific_interpretation
from app.interpretations.comparative_analysis import generate_comparative_summary

PURPLE = colors.HexColor("#7C4DFF")
DARK_BG = colors.HexColor("#0D1117")
DARK_CARD = colors.HexColor("#161B22")
BORDER_GRAY = colors.HexColor("#30363D")
TEXT_PRIMARY = colors.HexColor("#E8EAED")
TEXT_MUTED = colors.HexColor("#8B949E")
TEXT_DARK = colors.HexColor("#1A1A2E")
ACCENT_CYAN = colors.HexColor("#00BCD4")
ACCENT_EMERALD = colors.HexColor("#10B981")
ACCENT_AMBER = colors.HexColor("#F59E0B")
ACCENT_ROSE = colors.HexColor("#F43F5E")

NUC_COLORS_PDF = {
    "A": "#4CAF50",
    "U": "#F44336",
    "G": "#FF9800",
    "C": "#2196F3",
}


class CoverPage(Flowable):
    def __init__(self, analyses: list[AnalysisResult]):
        super().__init__()
        self.analyses = analyses
        self.width = letter[0]
        self.height = letter[1]

    def draw(self):
        c = self.canv
        w, h = self.width, self.height

        c.setFillColor(DARK_BG)
        c.rect(0, 0, w, h, fill=1, stroke=0)

        c.setFillColor(PURPLE, alpha=0.03)
        for i in range(5):
            cx = w * 0.5 + math.sin(i * 1.2) * 180
            cy = h * 0.6 + math.cos(i * 1.8) * 140
            c.circle(cx, cy, 120 + i * 30, fill=1, stroke=0)

        c.setStrokeColor(PURPLE, alpha=0.15)
        c.setLineWidth(0.5)
        for i in range(12):
            angle = i * math.pi / 6
            x1 = w / 2 + math.cos(angle) * 200
            y1 = h / 2 + math.sin(angle) * 200
            c.line(w / 2, h / 2, x1, y1)

        c.setFillColor(colors.white)
        c.setFont("Helvetica-Bold", 9)
        c.drawString(72, h - 46, "RNA GENOMICS INTELLIGENCE PLATFORM")

        c.setStrokeColor(PURPLE)
        c.setLineWidth(1.5)
        c.line(72, h - 50, 230, h - 50)

        c.setFillColor(PURPLE)
        c.setFont("Helvetica-Bold", 32)
        c.drawString(72, h - 110, "Therapeutic Candidate")

        c.setFillColor(colors.white)
        c.setFont("Helvetica-Bold", 32)
        c.drawString(72, h - 148, "Assessment Dossier")

        c.setFillColor(TEXT_MUTED)
        c.setFont("Helvetica", 11)
        c.drawString(72, h - 180, "Computational Thermodynamic & Structural Profiling Report")

        c.setStrokeColor(colors.HexColor("#1F2430"))
        c.setLineWidth(0.8)
        c.line(72, h - 198, w - 72, h - 198)

        meta_items = []
        if self.analyses:
            top = max(self.analyses, key=lambda a: _safe_val(a.stability.therapeutic_stability_score))
            meta_items.append(("Lead Candidate", top.sequence_name))
            meta_items.append(("Candidates Analyzed", str(len(self.analyses))))
            scores = [_safe_val(a.stability.therapeutic_stability_score) for a in self.analyses]
            meta_items.append(("Stability Range",
                               f"{min(scores):.1f} - {max(scores):.1f}"))
            if len(self.analyses) == 1:
                a0 = self.analyses[0]
                meta_items.append(("GC Content", f"{_safe_val(a0.stability.gc_content) * 100:.1f}%"))
                meta_items.append(("MFE", f"{_safe_val(a0.structure.minimum_free_energy):.2f} kcal/mol"))
                meta_items.append(("Structural Confidence", f"{_safe_val(a0.structure.structural_confidence) * 100:.1f}%"))
        meta_items.append(("Analysis Date", datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")))
        meta_items.append(("Report Reference", f"RNA-DOS-{datetime.utcnow().strftime('%Y%m%d-%H%M%S')}"))
        meta_items.append(("Document Classification", "Confidential — Research Use Only"))

        c.setFillColor(TEXT_MUTED)
        c.setFont("Helvetica", 9)
        y = h - 240
        for label, value in meta_items:
            c.setFillColor(TEXT_MUTED)
            c.drawString(72, y, label)
            c.setFillColor(colors.white)
            c.drawString(200, y, value)
            y -= 18

        c.setFillColor(TEXT_MUTED, alpha=0.4)
        c.setFont("Helvetica-Oblique", 8)
        c.drawCentredString(w / 2, 45, "CONFIDENTIALITY NOTICE")
        c.setFont("Helvetica", 7)
        c.drawCentredString(w / 2, 32, "This dossier contains proprietary computational genomics analyses intended solely for")
        c.drawCentredString(w / 2, 22, "therapeutic research and development purposes. Distribution beyond authorized personnel is prohibited.")

        c.setFont("Helvetica", 7)
        c.setFillColor(TEXT_MUTED, alpha=0.3)
        c.drawRightString(w - 72, 45, "Page 1")


class StructureTopologyFlowable(Flowable):
    """Enhanced RNA fold diagram that mirrors the SVG appearance."""

    def __init__(self, sequence, dot_bracket, base_pairs, width=440, height=260):
        super().__init__()
        self.sequence = sequence
        self.dot_bracket = dot_bracket
        self.base_pairs = base_pairs
        self.width = width
        self.height = height

    def draw(self):
        from app.utils.rna_structure_utils import parse_dot_bracket as _parse
        c = self.canv
        coords = generate_2d_coordinates(
            self.sequence, self.dot_bracket, self.base_pairs,
            width=int(self.width * 2), height=int(self.height * 2)
        )
        n = len(self.sequence)
        if n == 0 or not coords:
            return

        pair_map = {}
        for i, j in self.base_pairs:
            pair_map[i] = j
            pair_map[j] = i

        parsed = _parse(self.dot_bracket)

        xs = [p[0] for p in coords]
        ys = [p[1] for p in coords]
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        if max_x - min_x < 1:
            max_x = min_x + 100
        if max_y - min_y < 1:
            max_y = min_y + 100

        pad = 16
        plot_w = self.width - 2 * pad
        plot_h = self.height - 2 * pad
        target_fill = 0.82
        scale = min(
            plot_w * target_fill / (max_x - min_x),
            plot_h * target_fill / (max_y - min_y),
            3.5,
        )
        ox = pad + (plot_w - (max_x - min_x) * scale) / 2
        oy = pad + (plot_h - (max_y - min_y) * scale) / 2

        def tx(x):
            return ox + (x - min_x) * scale
        def ty(y):
            return oy + (y - min_y) * scale

        c.setFillColor(colors.HexColor("#0D1117"))
        c.roundRect(0, 0, self.width, self.height, 4, fill=1, stroke=0)
        c.setStrokeColor(colors.HexColor("#1F2430"))
        c.setLineWidth(0.5)
        c.roundRect(0, 0, self.width, self.height, 4, fill=0, stroke=1)

        # Stem highlights (subtle purple bars behind stems)
        for stem in parsed["stems"]:
            s_pairs = stem["pairs"]
            if s_pairs:
                s0 = s_pairs[0][0]
                s1 = s_pairs[-1][1]
                if s0 < len(coords) and s1 < len(coords):
                    x1 = tx(coords[s0][0]) - 4
                    x2 = tx(coords[s1][0]) + 4
                    y1 = min(ty(coords[s0][0]), ty(coords[s0][1])) - 3
                    y2 = max(ty(coords[s1][0]), ty(coords[s1][1])) + 3
                    if x2 > x1 and y2 > y1:
                        c.setFillColor(PURPLE, alpha=0.06)
                        c.roundRect(x1, y1, x2 - x1, y2 - y1, 2, fill=1, stroke=0)

        # Loop region indicators
        for loop in parsed["loops"]:
            ls, le = loop["start"], loop["end"]
            if ls < len(coords) and le < len(coords):
                lxs = [tx(coords[i][0]) for i in range(ls, min(le + 1, len(coords)))]
                lys = [ty(coords[i][1]) for i in range(ls, min(le + 1, len(coords)))]
                if lxs:
                    rx = (max(lxs) - min(lxs)) / 2 + 10
                    ry = (max(lys) - min(lys)) / 2 + 10
                    if rx > 0 and ry > 0:
                        c.setFillColor(ACCENT_CYAN, alpha=0.04)
                        c.setStrokeColor(ACCENT_CYAN, alpha=0.06)
                        c.setLineWidth(0.3)
                        c.setDash(3, 2)
                        cx_e = (min(lxs) + max(lxs)) / 2
                        cy_e = (min(lys) + max(lys)) / 2
                        p = c.beginPath()
                        p.ellipse(cx_e - rx, cy_e - ry, cx_e + rx, cy_e + ry)
                        c.drawPath(p, fill=1, stroke=1)
                        c.setDash()

        # Base-pair arcs
        for i, j in self.base_pairs:
            if i < j and i < len(coords) and j < len(coords):
                x1, y1 = tx(coords[i][0]), ty(coords[i][1])
                x2, y2 = tx(coords[j][0]), ty(coords[j][1])
                dx = x2 - x1
                dy = y2 - y1
                dist = math.sqrt(dx * dx + dy * dy)
                if dist > 3:
                    mid_x = (x1 + x2) / 2
                    mid_y = (y1 + y2) / 2
                    perp_x = -dy / dist * 16
                    perp_y = dx / dist * 16
                    cp_x = mid_x + perp_x
                    cp_y = mid_y + perp_y - abs(perp_y) * 0.3
                    c.setStrokeColor(PURPLE, alpha=0.25)
                    c.setLineWidth(0.4)
                    c.setDash(2, 2)
                    p = c.beginPath()
                    p.moveTo(x1, y1)
                    p.curveTo(cp_x, cp_y, cp_x, cp_y, x2, y2)
                    c.drawPath(p, fill=0, stroke=1)
                    c.setDash()

        # Backbone
        c.setStrokeColor(colors.HexColor("#30363D"))
        c.setLineWidth(0.8)
        for i in range(n - 1):
            if i < len(coords) and i + 1 < len(coords):
                x1, y1 = tx(coords[i][0]), ty(coords[i][1])
                x2, y2 = tx(coords[i + 1][0]), ty(coords[i + 1][1])
                c.line(x1, y1, x2, y2)

        # Nucleotide nodes
        for i in range(n):
            if i >= len(coords):
                break
            x, y = tx(coords[i][0]), ty(coords[i][1])
            nt = self.sequence[i] if i < len(self.sequence) else "?"
            color_hex = NUC_COLORS_PDF.get(nt, "#666")
            clr = colors.HexColor(color_hex)
            is_paired = i in pair_map
            r = 5 if is_paired else 3.8

            if is_paired:
                c.setFillColor(clr, alpha=0.12)
                c.circle(x, y, r + 3.5, fill=1, stroke=0)

            c.setFillColor(clr)
            c.setStrokeColor(colors.white, alpha=0.3 if is_paired else 0.08)
            c.setLineWidth(0.3)
            c.circle(x, y, r, fill=1, stroke=1)

            font_size = 5.5 if is_paired else 4.5
            c.setFont("Helvetica-Bold", font_size)
            c.setFillColor(colors.white)
            c.drawCentredString(x, y - 2, nt)

        c.setFont("Helvetica-Bold", 7.5)
        c.setFillColor(colors.HexColor("#C9D1D9"))
        c.drawCentredString(self.width / 2, self.height - 7,
                            "RNA Secondary Structure Topology — True Fold Layout")

        num_pairs = len(self.base_pairs)
        paired_pct = (2 * num_pairs / n * 100) if n > 0 else 0
        stem_count = len(parsed["stems"])
        c.setFont("Helvetica", 5.5)
        c.setFillColor(TEXT_MUTED)
        c.drawRightString(self.width - 6, 3.5,
                          f"{n} nt | {num_pairs} bp | {stem_count} stems | Paired: {paired_pct:.1f}%")


class GaugeBar(Flowable):
    def __init__(self, label, value, max_value=100.0, width=440, height=22, color=PURPLE):
        super().__init__()
        self.label = label
        self.value = value
        self.max_value = max_value
        self.width = width
        self.height = height
        self.color = color

    def draw(self):
        c = self.canv
        frac = min(1.0, self.value / self.max_value)
        bar_x = 120
        bar_w = self.width - bar_x
        bar_h = 12
        bar_y = (self.height - bar_h) / 2

        c.setFillColor(colors.HexColor("#1A1F2E"))
        c.roundRect(bar_x, bar_y, bar_w, bar_h, 3, fill=1, stroke=0)

        fill_color = self.color
        if frac < 0.33:
            fill_color = ACCENT_ROSE
        elif frac < 0.66:
            fill_color = ACCENT_AMBER

        c.setFillColor(fill_color)
        if frac > 0.01:
            c.roundRect(bar_x, bar_y, bar_w * frac, bar_h, 3, fill=1, stroke=0)

        c.setStrokeColor(colors.HexColor("#1F2430"))
        c.setLineWidth(0.5)
        c.roundRect(bar_x, bar_y, bar_w, bar_h, 3, fill=0, stroke=1)

        c.setFont("Helvetica-Bold", 7.5)
        c.setFillColor(colors.HexColor("#C9D1D9"))
        c.drawString(0, bar_y + 1, self.label)

        c.setFont("Helvetica-Bold", 8.5)
        c.setFillColor(colors.white)
        c.drawRightString(self.width, bar_y + 1, f"{self.value:.1f}")


class GaugeMeter(Flowable):
    def __init__(self, label, value, unit="", width=200, height=70, color=PURPLE):
        super().__init__()
        self.label = label
        self.value = value
        self.unit = unit
        self.width = width
        self.height = height
        self.color = color

    def draw(self):
        c = self.canv
        cx = self.width / 2
        r = 26

        c.setFillColor(colors.HexColor("#0D1117"))
        c.roundRect(0, 0, self.width, self.height, 4, fill=1, stroke=0)
        c.setStrokeColor(colors.HexColor("#1F2430"))
        c.setLineWidth(0.5)
        c.roundRect(0, 0, self.width, self.height, 4, fill=0, stroke=1)

        c.setStrokeColor(colors.HexColor("#1F2430"))
        c.setLineWidth(3)
        c.circle(cx, self.height - 30, r, fill=0, stroke=1)

        frac = min(1.0, max(0.0, self.value / 100.0))
        angle = -180 * frac
        start_angle = -180
        end_angle = start_angle + angle

        c.setStrokeColor(self.color)
        c.setLineWidth(3)
        c.arc(cx - r, self.height - 30 - r,
              cx + r, self.height - 30 + r,
              start_angle, end_angle)

        c.setFont("Helvetica-Bold", 11)
        c.setFillColor(colors.white)
        c.drawCentredString(cx, self.height - 34, f"{self.value:.0f}")

        if self.unit:
            c.setFont("Helvetica", 6)
            c.setFillColor(TEXT_MUTED)
            c.drawCentredString(cx, self.height - 46, self.unit)

        c.setFont("Helvetica", 6.5)
        c.setFillColor(colors.HexColor("#8B949E"))
        c.drawCentredString(cx, 8, self.label)


class CompBarChart(Flowable):
    def __init__(self, data, title="", width=440, height=160):
        super().__init__()
        self.data = data
        self.title = title
        self.width = width
        self.height = height

    def draw(self):
        c = self.canv
        items = list(self.data.items())
        n = len(items)
        if n == 0:
            return

        margin_l = 50
        margin_r = 20
        margin_t = 18
        margin_b = 30
        plot_w = self.width - margin_l - margin_r
        plot_h = self.height - margin_t - margin_b
        bar_w = min(30, plot_w / n * 0.6)
        gap = plot_w / n

        c.setFillColor(colors.HexColor("#0D1117"))
        c.roundRect(0, 0, self.width, self.height, 4, fill=1, stroke=0)
        c.setStrokeColor(colors.HexColor("#1F2430"))
        c.setLineWidth(0.5)
        c.roundRect(0, 0, self.width, self.height, 4, fill=0, stroke=1)

        c.setFont("Helvetica-Bold", 8)
        c.setFillColor(colors.HexColor("#C9D1D9"))
        c.drawCentredString(self.width / 2, self.height - 4, self.title)

        max_val = max(v for _, v in items) if items else 1

        for i, (label, val) in enumerate(items):
            x = margin_l + i * gap + (gap - bar_w) / 2
            bar_h = (val / max_val) * plot_h if max_val > 0 else 0
            y = margin_b

            colors_map = {
                "A": colors.HexColor("#4CAF50"),
                "C": colors.HexColor("#2196F3"),
                "G": colors.HexColor("#FF9800"),
                "U": colors.HexColor("#F44336"),
            }
            clr = colors_map.get(label, PURPLE)
            c.setFillColor(clr)
            c.rect(x, y, bar_w, bar_h, fill=1, stroke=0)

            c.setFont("Helvetica-Bold", 7)
            c.setFillColor(colors.white)
            c.drawCentredString(x + bar_w / 2, y + bar_h + 2, f"{val:.1f}%")

            c.setFont("Helvetica", 7)
            c.setFillColor(TEXT_MUTED)
            c.drawCentredString(x + bar_w / 2, 6, label)


class RadarChart(Flowable):
    def __init__(self, categories, values, title="", width=240, height=240):
        super().__init__()
        self.categories = categories
        self.values = values
        self.title = title
        self.width = width
        self.height = height

    def draw(self):
        c = self.canv
        n = len(self.categories)
        if n < 3:
            return

        cx = self.width / 2
        cy = self.height / 2 - 10
        r = min(cx, cy) - 20

        c.setFillColor(colors.HexColor("#0D1117"))
        c.roundRect(0, 0, self.width, self.height, 4, fill=1, stroke=0)
        c.setStrokeColor(colors.HexColor("#1F2430"))
        c.setLineWidth(0.5)
        c.roundRect(0, 0, self.width, self.height, 4, fill=0, stroke=1)

        for level in [0.25, 0.5, 0.75, 1.0]:
            pts = []
            for i in range(n):
                angle = -math.pi / 2 + 2 * math.pi * i / n
                pts.append((cx + r * level * math.cos(angle),
                            cy + r * level * math.sin(angle)))
            c.setStrokeColor(colors.HexColor("#1F2430"))
            c.setLineWidth(0.5)
            p = c.beginPath()
            p.moveTo(*pts[0])
            for pt in pts[1:]:
                p.lineTo(*pt)
            p.close()
            c.drawPath(p, fill=0, stroke=1)

        for i in range(n):
            angle = -math.pi / 2 + 2 * math.pi * i / n
            c.line(cx, cy, cx + r * math.cos(angle), cy + r * math.sin(angle))

        val_pts = []
        for i in range(n):
            v = max(0, min(1, self.values[i] / 100.0 if n else 0))
            angle = -math.pi / 2 + 2 * math.pi * i / n
            val_pts.append((cx + r * v * math.cos(angle),
                            cy + r * v * math.sin(angle)))

        if val_pts:
            c.setFillColor(PURPLE, alpha=0.12)
            c.setStrokeColor(PURPLE)
            c.setLineWidth(1.5)
            p = c.beginPath()
            p.moveTo(*val_pts[0])
            for pt in val_pts[1:]:
                p.lineTo(*pt)
            p.close()
            c.drawPath(p, fill=1, stroke=1)

            for pt in val_pts:
                c.setFillColor(PURPLE)
                c.circle(pt[0], pt[1], 3, fill=1, stroke=0)

        for i in range(n):
            angle = -math.pi / 2 + 2 * math.pi * i / n
            label_x = cx + (r + 18) * math.cos(angle)
            label_y = cy + (r + 18) * math.sin(angle)
            c.setFont("Helvetica", 6)
            c.setFillColor(colors.HexColor("#8B949E"))
            c.drawCentredString(label_x, label_y - 2.5, self.categories[i])

        c.setFont("Helvetica-Bold", 7.5)
        c.setFillColor(colors.HexColor("#C9D1D9"))
        c.drawCentredString(self.width / 2, self.height - 4, self.title)


class RiskBadge(Flowable):
    def __init__(self, risk_type, level, description, width=440):
        super().__init__()
        self.risk_type = risk_type
        self.level = level
        self.description = description
        self.width = width
        self.height = 28

    def draw(self):
        c = self.canv
        clr_map = {"low": ACCENT_EMERALD, "moderate": ACCENT_AMBER, "high": ACCENT_ROSE, "critical": colors.HexColor("#FF1744")}
        level = self.level.lower()
        clr = clr_map.get(level, ACCENT_AMBER)

        c.setFillColor(colors.HexColor("#0D1117"))
        c.roundRect(0, 0, self.width, self.height, 3, fill=1, stroke=0)
        c.setStrokeColor(clr, alpha=0.3)
        c.setLineWidth(0.5)
        c.roundRect(0, 0, self.width, self.height, 3, fill=0, stroke=1)

        c.setFillColor(clr)
        c.circle(14, self.height / 2, 4.5, fill=1, stroke=0)

        c.setFont("Helvetica-Bold", 7.5)
        c.setFillColor(clr)
        c.drawString(24, self.height / 2 - 3, self.risk_type.upper())

        c.setFont("Helvetica", 7)
        c.setFillColor(colors.HexColor("#C9D1D9"))
        c.drawString(90, self.height / 2 - 3, self.description)

        c.setFont("Helvetica-Bold", 7)
        c.setFillColor(clr)
        label_w = c.stringWidth(self.level.upper(), "Helvetica-Bold", 7)
        c.drawRightString(self.width - 10, self.height / 2 - 3, self.level.upper())


class StabilityRadarFlowable(Flowable):
    def __init__(self, analyses, width=440, height=200):
        super().__init__()
        self.analyses = analyses
        self.width = width
        self.height = height

    def draw(self):
        c = self.canv
        n_analyses = len(self.analyses)
        if n_analyses == 0:
            return

        margin_l = 50
        margin_r = 20
        margin_t = 20
        margin_b = 30
        plot_w = self.width - margin_l - margin_r
        plot_h = self.height - margin_t - margin_b

        c.setFillColor(colors.HexColor("#0D1117"))
        c.roundRect(0, 0, self.width, self.height, 4, fill=1, stroke=0)
        c.setStrokeColor(colors.HexColor("#1F2430"))
        c.setLineWidth(0.5)
        c.roundRect(0, 0, self.width, self.height, 4, fill=0, stroke=1)

        c.setFont("Helvetica-Bold", 7.5)
        c.setFillColor(colors.HexColor("#C9D1D9"))
        c.drawCentredString(self.width / 2, self.height - 4,
                            "Candidate Stability Comparison")

        metrics = [
            ("Stability", lambda a: a.stability.therapeutic_stability_score, PURPLE),
            ("GC%", lambda a: a.stability.gc_content * 100, ACCENT_CYAN),
            ("Confidence", lambda a: a.structure.structural_confidence * 100, ACCENT_EMERALD),
            ("Pairing", lambda a: a.structure.paired_fraction * 100, ACCENT_AMBER),
        ]

        n_metrics = len(metrics)
        bar_w = min(20, plot_w / (n_metrics * n_analyses) * 0.8)
        group_gap = plot_w / n_metrics

        max_val = 100

        for mi, (mname, mfunc, mcolor) in enumerate(metrics):
            group_x = margin_l + mi * group_gap
            for ai, a in enumerate(self.analyses):
                val = mfunc(a)
                bar_h = (val / max_val) * plot_h
                x = group_x + (group_gap - bar_w * n_analyses) / 2 + ai * bar_w
                y = margin_b

                colors_list = [PURPLE, ACCENT_CYAN, ACCENT_EMERALD, ACCENT_AMBER, ACCENT_ROSE]
                clr = colors_list[ai % len(colors_list)]
                c.setFillColor(clr, alpha=0.7)
                c.rect(x, y, bar_w, bar_h, fill=1, stroke=0)

            c.setFont("Helvetica", 6)
            c.setFillColor(TEXT_MUTED)
            c.drawCentredString(group_x + group_gap / 2, 6, mname)

        c.setFont("Helvetica", 6)
        c.setFillColor(TEXT_MUTED)
        for ai, a in enumerate(self.analyses):
            colors_list = [PURPLE, ACCENT_CYAN, ACCENT_EMERALD, ACCENT_AMBER, ACCENT_ROSE]
            clr = colors_list[ai % len(colors_list)]
            x = margin_l + plot_w - 100 + ai * 50
            c.setFillColor(clr, alpha=0.7)
            c.rect(x, self.height - 14, 8, 8, fill=1, stroke=0)
            c.setFillColor(TEXT_MUTED)
            c.setFont("Helvetica", 5.5)
            c.drawString(x + 10, self.height - 12, a.sequence_name[:12])


def _safe_val(val, default=0.0):
    """Return a safe numeric value, replacing NaN/None with default."""
    if val is None:
        return default
    try:
        return val if math.isfinite(val) else default
    except (TypeError, ValueError):
        return default


def _generate_executive_summary(a):
    gc_pct = _safe_val(a.stability.gc_content) * 100
    score = _safe_val(a.stability.therapeutic_stability_score)
    mfe = _safe_val(a.structure.minimum_free_energy)
    conf = _safe_val(a.structure.structural_confidence) * 100
    pairs = len(a.structure.base_pairs or [])
    paired_pct = _safe_val(a.structure.paired_fraction) * 100
    stem_metrics = compute_stem_metrics(a.structure.dot_bracket or "", a.structure.base_pairs or [])

    if score >= 70:
        verdict = "strong therapeutic candidate"
        stability_note = (
            "Robust thermodynamic persistence with favorable folding "
            "energetics supporting conformational integrity under "
            "physiological conditions."
        )
    elif score >= 45:
        verdict = "moderate therapeutic candidate"
        stability_note = (
            "Acceptable structural stability; sequence optimization "
            "may enhance conformational rigidity and nuclease resistance."
        )
    else:
        verdict = "limited therapeutic candidate"
        stability_note = (
            "Marginal structural persistence; significant engineering "
            "required for therapeutic-grade conformational stability."
        )

    summary = (
        f"Candidate {a.sequence_name} presents as a {verdict} "
        f"(TSI: {score:.1f}/100). "
        f"{stability_note} "
        f"Fold comprises {pairs} canonical base pairs "
        f"({paired_pct:.1f}% paired) at MFE {mfe:.2f} kcal/mol. "
        f"Structural confidence: {conf:.1f}%, GC content: {gc_pct:.1f}%. "
    )

    if stem_metrics["stem_count"] > 0:
        summary += (
            f"Topology: {stem_metrics['stem_count']} stem(s), "
            f"avg {stem_metrics['avg_stem_length']:.1f} bp, "
            f"{stem_metrics['loop_count']} loop(s), "
            f"avg {stem_metrics['avg_loop_size']:.1f} nt. "
        )

    if a.stability.ensemble_diversity < 0.2:
        summary += (
            "Low ensemble diversity confirms a well-defined structural "
            "scaffold with minimal conformational heterogeneity."
        )
    elif a.stability.ensemble_diversity > 0.5:
        summary += (
            "Elevated ensemble diversity indicates conformational "
            "plasticity that may influence target-binding kinetics."
        )
    else:
        summary += (
            "Moderate ensemble diversity reflects a balanced "
            "conformational landscape with a dominant fold population."
        )

    return summary, verdict


def _generate_thermo_interpretation(a):
    mfe = _safe_val(a.structure.minimum_free_energy)
    mfe_per_nt = _safe_val(a.stability.mfe_per_nucleotide)
    diversity = _safe_val(a.stability.ensemble_diversity)
    paired_pct = _safe_val(a.structure.paired_fraction) * 100
    gc_pct = _safe_val(a.stability.gc_content) * 100
    score = _safe_val(a.stability.therapeutic_stability_score)

    parts = []

    parts.append(f"MFE: {mfe:.2f} kcal/mol | Normalized: {mfe_per_nt:.2f} kcal/mol·nt")
    parts.append(
        f"Thermodynamic favorability: "
        + (
            f"strong base-pair stacking, well-defined fold."
            if mfe_per_nt <= -0.3
            else f"moderate stability; competing low-energy folds possible."
            if mfe_per_nt <= -0.15
            else f"shallow landscape; multiple near-ground-state conformations."
        )
    )

    parts.append(f"GC Content: {gc_pct:.1f}%")
    parts.append(
        (
            f"Enhanced H-bond density (+3 per G-C) and Pi-stacking "
            f"elevate Tm and nuclease resistance."
            if gc_pct >= 55
            else f"Balanced stabilization with conformational flexibility "
                 f"for target engagement."
            if gc_pct >= 40
            else f"Reduced H-bond density; lower thermal stability, "
                 f"more dynamic conformational behavior."
        )
    )

    parts.append(f"Ensemble Diversity: {diversity:.3f}")
    parts.append(
        (
            f"Dominant fold > 80% population; high conformational homogeneity."
            if diversity < 0.2
            else f"Primary fold with subpopulations; minor conformational switching."
            if diversity < 0.5
            else f"Highly heterogeneous ensemble; multiple nearly-isoenergetic states."
        )
    )

    parts.append(f"Paired Fraction: {paired_pct:.1f}%")
    parts.append(
        (
            f"Compact fold; may occlude target-recognition motifs."
            if paired_pct >= 60
            else f"Ordered structure with accessible single-stranded regions."
            if paired_pct >= 35
            else f"Predominantly single-stranded; limited structural ordering."
        )
    )

    parts.append(f"Therapeutic Stability Score: {score:.1f}/100")
    parts.append(
        (
            f"Excellent profile — candidate suitable for preclinical advancement."
            if score >= 70
            else f"Acceptable — targeted sequence optimization recommended."
            if score >= 45
            else f"Suboptimal — iterative redesign and chemical modification required."
        )
    )

    return "\n\n".join(parts)


def _generate_scientific_interpretation(a):
    parts = []
    score = _safe_val(a.stability.therapeutic_stability_score)
    gc_pct = _safe_val(a.stability.gc_content) * 100
    mfe = _safe_val(a.structure.minimum_free_energy)
    diversity = _safe_val(a.stability.ensemble_diversity)
    paired_pct = _safe_val(a.structure.paired_fraction) * 100
    stem_metrics = compute_stem_metrics(a.structure.dot_bracket or "", a.structure.base_pairs or [])
    conf = _safe_val(a.structure.structural_confidence) * 100

    parts.append("STRUCTURAL STABILITY")
    parts.append(
        f"Therapeutic Stability Index {score:.1f}/100 — "
        f"MFE {mfe:.2f} kcal/mol, base-pairing density "
        f"{paired_pct:.1f}%, ensemble diversity {diversity:.3f}. "
        + (
            f"Thermodynamically favored fold with well-defined energy landscape."
            if score >= 70
            else f"Metastable conformation; sequence optimization may improve rigidity."
            if score >= 45
            else f"Marginal stability; substantial sequence engineering required."
        )
    )

    parts.append("DEGRADATION RESISTANCE")
    parts.append(
        f"GC content {gc_pct:.1f}% — "
        + (
            f"elevated G-C triplet hydrogen bonding enhances nuclease resistance and Tm."
            if gc_pct >= 55
            else f"moderate stabilization; modified nucleotides (2'-OMe, PS) may augment half-life."
            if gc_pct >= 40
            else f"AU-rich composition predisposes to endonuclease cleavage. Chemical modification recommended."
        )
    )

    parts.append("STRUCTURAL COMPACTNESS")
    parts.append(
        f"Pairing density {paired_pct:.1f}% across "
        f"{stem_metrics['stem_count']} stem(s), {stem_metrics['loop_count']} loop(s). "
        + (
            f"Extensive secondary structure; verify start codon accessibility."
            if paired_pct >= 60
            else f"Balanced ordering with accessible single-stranded regions."
            if paired_pct >= 35
            else f"Open conformation; limited structural persistence."
        )
    )

    parts.append("THERAPEUTIC VIABILITY")
    parts.append(
        f"Candidate classified as "
        + (
            f"Tier I (high priority). Ensemble homogeneity supports predictable in vivo behavior."
            if score >= 70
            else f"Tier II (viable). Moderate optimization recommended."
            if score >= 45
            else f"Tier III (suboptimal). Sequence redesign and chemical stabilization required."
        )
    )

    parts.append("DELIVERY PROFILE")
    length = _safe_val(a.length)
    if length > 100:
        parts.append(f"Length {length:.0f} nt exceeds siRNA range; LNP or GalNAc delivery recommended.")
    elif length > 21:
        parts.append(f"Length {length:.0f} nt compatible with LNP and GalNAc platforms.")
    else:
        parts.append(f"Short sequence ({length:.0f} nt) amenable to standard transfection.")
    if gc_pct > 65:
        parts.append("High GC content may affect solubility; pre-formulation characterization advised.")

    return "\n\n".join(parts)


def _generate_risk_assessment(a):
    risks = []
    gc_pct = _safe_val(a.stability.gc_content) * 100
    paired_pct = _safe_val(a.structure.paired_fraction) * 100
    conf = _safe_val(a.structure.structural_confidence)
    diversity = _safe_val(a.stability.ensemble_diversity)
    mfe_per_nt = _safe_val(a.stability.mfe_per_nucleotide)
    stem_metrics = compute_stem_metrics(a.structure.dot_bracket or "", a.structure.base_pairs or [])

    if gc_pct < 35:
        risks.append(("GC Imbalance", "high",
                      f"Low GC content ({gc_pct:.1f}%) may compromise thermal stability and nuclease resistance."))
    elif gc_pct < 45:
        risks.append(("GC Imbalance", "moderate",
                      f"Moderate GC content ({gc_pct:.1f}%) is acceptable but suboptimal for maximum stability."))
    else:
        risks.append(("GC Imbalance", "low",
                      f"GC content ({gc_pct:.1f}%) within optimal range for structural stability."))

    if paired_pct < 30:
        risks.append(("Low Pairing", "high",
                      f"Only {paired_pct:.1f}% of nucleotides are paired, indicating a predominantly unstructured RNA."))
    elif paired_pct < 50:
        risks.append(("Low Pairing", "moderate",
                      f"Pairing fraction ({paired_pct:.1f}%) suggests moderate structural organization."))
    else:
        risks.append(("Low Pairing", "low",
                      f"Pairing fraction ({paired_pct:.1f}%) indicates well-structured topology."))

    if mfe_per_nt > -0.15:
        risks.append(("Instability Indicator", "high",
                      f"Shallow MFE per nucleotide ({mfe_per_nt:.3f} kcal/mol/nt) suggests limited thermodynamic driving force for folding."))
    elif mfe_per_nt > -0.3:
        risks.append(("Instability Indicator", "moderate",
                      f"MFE per nucleotide ({mfe_per_nt:.3f} kcal/mol/nt) indicates moderate folding stability."))
    else:
        risks.append(("Instability Indicator", "low",
                      f"Favorable MFE per nucleotide ({mfe_per_nt:.3f} kcal/mol/nt) supports stable fold formation."))

    if conf < 0.5:
        risks.append(("Folding Irregularity", "high",
                      f"Low structural confidence ({conf * 100:.1f}%) suggests significant conformational heterogeneity."))
    elif conf < 0.75:
        risks.append(("Folding Irregularity", "moderate",
                      f"Moderate structural confidence ({conf * 100:.1f}%) indicates some conformational flexibility."))
    else:
        risks.append(("Folding Irregularity", "low",
                      f"High structural confidence ({conf * 100:.1f}%) supports a deterministic folding pathway."))

    if diversity > 0.5:
        risks.append(("Ensemble Diversity", "high",
                      f"Elevated ensemble diversity ({diversity:.3f}) indicates competing conformational states."))
    elif diversity > 0.25:
        risks.append(("Ensemble Diversity", "moderate",
                      f"Moderate ensemble diversity ({diversity:.3f}) suggests manageable conformational heterogeneity."))
    else:
        risks.append(("Ensemble Diversity", "low",
                      f"Low ensemble diversity ({diversity:.3f}) supports a single dominant fold."))

    if stem_metrics["stem_count"] == 0:
        risks.append(("Structural Organization", "high",
                      "No stem regions detected; the RNA may lack defined secondary structure."))
    elif stem_metrics["stem_count"] < 2:
        risks.append(("Structural Organization", "moderate",
                      f"Only {stem_metrics['stem_count']} stem region(s) detected; limited structural hierarchy."))
    else:
        risks.append(("Structural Organization", "low",
                      f"{stem_metrics['stem_count']} stem region(s) provide a well-organized structural scaffold."))

    return risks


def _build_composition_table(a):
    styles = getSampleStyleSheet()
    s_table_body = ParagraphStyle("tb", parent=styles["Normal"], fontName="Helvetica", fontSize=7.5, leading=10, textColor=colors.HexColor("#C9D1D9"))

    comp_headers = [
        Paragraph("<b>Nucleotide</b>", s_table_body),
        Paragraph("<b>Count</b>", s_table_body),
        Paragraph("<b>Percentage</b>", s_table_body),
        Paragraph("<b>Distribution</b>", s_table_body),
    ]
    seq_counts = {
        "A": a.sequence.count("A"),
        "C": a.sequence.count("C"),
        "G": a.sequence.count("G"),
        "U": a.sequence.count("U"),
    }
    total = len(a.sequence)
    nt_info = [
        ("A", "Adenine", colors.HexColor("#4CAF50")),
        ("C", "Cytosine", colors.HexColor("#2196F3")),
        ("G", "Guanine", colors.HexColor("#FF9800")),
        ("U", "Uracil", colors.HexColor("#F44336")),
    ]

    rows = [comp_headers]
    for key, label, clr in nt_info:
        cnt = seq_counts.get(key, 0)
        pct = (cnt / total * 100) if total > 0 else 0
        bar_len = max(1, int(pct / 2.5))
        bar = "█" * bar_len
        rows.append([
            Paragraph(f"<b>{key}</b> — {label}", s_table_body),
            Paragraph(str(cnt), s_table_body),
            Paragraph(f"{pct:.1f}%", s_table_body),
            Paragraph(f'<font color="{clr.hexval()}">{bar}</font>', s_table_body),
        ])

    t = Table(rows, colWidths=[130, 50, 70, 170])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), PURPLE),
        ("ALIGN", (0, 0), (-1, -1), "LEFT"),
        ("ALIGN", (1, 0), (2, -1), "CENTER"),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
        ("LINEBELOW", (0, 0), (-1, -1), 0.3, BORDER_GRAY),
        ("BACKGROUND", (0, 1), (-1, -1), colors.HexColor("#0D1117")),
    ]))
    return t


def _build_metrics_table(a):
    styles = getSampleStyleSheet()
    s_td = ParagraphStyle("td", parent=styles["Normal"], fontName="Helvetica", fontSize=7.5, leading=10, textColor=colors.HexColor("#C9D1D9"))
    s_td_val = ParagraphStyle("tdv", parent=styles["Normal"], fontName="Helvetica-Bold", fontSize=7.5, leading=10, textColor=colors.white)

    stem_metrics = compute_stem_metrics(a.structure.dot_bracket or "", a.structure.base_pairs or [])

    metrics_data = [
        [Paragraph("<b>Metric</b>", s_td), Paragraph("<b>Value</b>", s_td), Paragraph("<b>Metric</b>", s_td), Paragraph("<b>Value</b>", s_td)],
        [Paragraph("Therapeutic Score", s_td), Paragraph(f"{a.stability.therapeutic_stability_score:.2f} / 100", s_td_val),
         Paragraph("GC Content", s_td), Paragraph(f"{a.stability.gc_content * 100:.1f}%", s_td_val)],
        [Paragraph("Minimum Free Energy", s_td), Paragraph(f"{a.structure.minimum_free_energy:.2f} kcal/mol", s_td_val),
         Paragraph("MFE / Nucleotide", s_td), Paragraph(f"{a.stability.mfe_per_nucleotide:.4f} kcal/mol/nt", s_td_val)],
        [Paragraph("Structural Confidence", s_td), Paragraph(f"{a.structure.structural_confidence * 100:.1f}%", s_td_val),
         Paragraph("Ensemble Diversity", s_td), Paragraph(f"{a.stability.ensemble_diversity:.4f}", s_td_val)],
        [Paragraph("Paired Fraction", s_td), Paragraph(f"{a.structure.paired_fraction * 100:.1f}%", s_td_val),
         Paragraph("Base Pair Count", s_td), Paragraph(str(len(a.structure.base_pairs)), s_td_val)],
        [Paragraph("Stem Regions", s_td), Paragraph(str(stem_metrics["stem_count"]), s_td_val),
         Paragraph("Avg Stem Length", s_td), Paragraph(f"{stem_metrics['avg_stem_length']:.1f} bp", s_td_val)],
        [Paragraph("Loop Regions", s_td), Paragraph(str(stem_metrics["loop_count"]), s_td_val),
         Paragraph("Avg Loop Size", s_td), Paragraph(f"{stem_metrics['avg_loop_size']:.1f} nt", s_td_val)],
        [Paragraph("Sequence Length", s_td), Paragraph(f"{a.length} nt", s_td_val),
         Paragraph("Max Pairing Depth", s_td), Paragraph(str(stem_metrics["max_depth"]), s_td_val)],
    ]

    t = Table(metrics_data, colWidths=[100, 120, 100, 120])
    t.setStyle(TableStyle([
        ("ALIGN", (0, 0), (-1, -1), "LEFT"),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("GRID", (0, 0), (-1, -1), 0.3, BORDER_GRAY),
        ("BACKGROUND", (0, 0), (-1, 0), PURPLE),
        ("BACKGROUND", (1, 1), (1, -1), colors.HexColor("#0A0E17")),
        ("BACKGROUND", (3, 1), (3, -1), colors.HexColor("#0A0E17")),
        ("BACKGROUND", (0, 1), (0, -1), colors.HexColor("#0D1117")),
        ("BACKGROUND", (2, 1), (2, -1), colors.HexColor("#0D1117")),
    ]))
    return t


def generate_pdf_report(analyses: list[AnalysisResult]) -> bytes:
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(
        buffer,
        pagesize=letter,
        leftMargin=54,
        rightMargin=54,
        topMargin=60,
        bottomMargin=60,
    )

    styles = getSampleStyleSheet()

    s_section = ParagraphStyle(
        "SectionHeading", parent=styles["Normal"],
        fontName="Helvetica-Bold", fontSize=14, leading=18,
        textColor=PURPLE, spaceBefore=14, spaceAfter=8, keepWithNext=True,
    )
    s_subsection = ParagraphStyle(
        "SubHeading", parent=styles["Normal"],
        fontName="Helvetica-Bold", fontSize=10.5, leading=14,
        textColor=colors.HexColor("#E8EAED"), spaceBefore=12, spaceAfter=6,
    )
    s_subsubsection = ParagraphStyle(
        "SubSubHeading", parent=styles["Normal"],
        fontName="Helvetica-Bold", fontSize=9, leading=12,
        textColor=ACCENT_CYAN, spaceBefore=10, spaceAfter=4,
    )
    s_body = ParagraphStyle(
        "ReportBody", parent=styles["Normal"],
        fontName="Helvetica", fontSize=8.5, leading=12,
        textColor=colors.HexColor("#C9D1D9"), spaceAfter=6,
    )
    s_body_small = ParagraphStyle(
        "ReportBodySmall", parent=styles["Normal"],
        fontName="Helvetica", fontSize=7.5, leading=10,
        textColor=colors.HexColor("#C9D1D9"), spaceAfter=4,
    )
    s_code = ParagraphStyle(
        "CodeStyle", parent=styles["Normal"],
        fontName="Courier", fontSize=7, leading=9,
        textColor=ACCENT_CYAN, spaceAfter=4,
    )
    s_caption = ParagraphStyle(
        "Caption", parent=styles["Normal"],
        fontName="Helvetica-Oblique", fontSize=7, leading=9,
        textColor=TEXT_MUTED, spaceAfter=4,
    )
    s_table_header = ParagraphStyle(
        "TableHeader", parent=styles["Normal"],
        fontName="Helvetica-Bold", fontSize=7, leading=9,
        textColor=colors.white,
    )
    s_table_body = ParagraphStyle(
        "TableBody", parent=styles["Normal"],
        fontName="Helvetica", fontSize=7, leading=9,
        textColor=colors.HexColor("#C9D1D9"),
    )

    story: list[Any] = []

    # ================================================================== #
    # COVER PAGE
    # ================================================================== #
    story.append(CoverPage(analyses))
    story.append(PageBreak())

    # ================================================================== #
    # SECTION 1 — EXECUTIVE SUMMARY
    # ================================================================== #
    story.append(Paragraph("1. Executive Summary", s_section))
    story.append(HRFlowable(width="100%", thickness=0.5, color=PURPLE, spaceAfter=8))

    if len(analyses) == 1:
        a = analyses[0]
        summary_text, verdict = _generate_executive_summary(a)
        story.append(Paragraph(summary_text, s_body))
        story.append(Spacer(1, 4))

        rec_color = ACCENT_EMERALD if a.stability.therapeutic_stability_score >= 60 else ACCENT_AMBER
        rec_text = (
            "Proceed to preclinical evaluation and sequence optimization studies."
            if a.stability.therapeutic_stability_score >= 60
            else "Consider sequence engineering to enhance structural stability before therapeutic assessment."
        )
        story.append(Paragraph(
            f'<b><font color="{rec_color.hexval()}">Recommendation:</font></b> {rec_text}',
            s_body_small,
        ))
    else:
        avg_score = sum(a.stability.therapeutic_stability_score for a in analyses) / len(analyses)
        top = max(analyses, key=lambda x: x.stability.therapeutic_stability_score)
        story.append(Paragraph(
            f"This dossier presents a comparative assessment of {len(analyses)} "
            f"RNA candidate sequences. The cohort exhibits a mean Therapeutic "
            f"Stability Index of {avg_score:.1f}/100, with candidate "
            f"'{top.sequence_name}' achieving the highest score "
            f"({top.stability.therapeutic_stability_score:.1f}/100). "
            f"Candidate ranking is based on integrated thermodynamic stability, "
            f"structural confidence, and compositional analysis.",
            s_body,
        ))

    story.append(Spacer(1, 12))

    # ================================================================== #
    # SECTION 2 — THERMODYNAMIC ANALYSIS (summary across all candidates)
    # ================================================================== #
    story.append(Paragraph("2. Thermodynamic &amp; Structural Overview", s_section))
    story.append(HRFlowable(width="100%", thickness=0.5, color=PURPLE, spaceAfter=8))

    if len(analyses) >= 2:
        story.append(StabilityRadarFlowable(analyses, width=460, height=180))
        story.append(Spacer(1, 4))
        story.append(Paragraph(
            "Figure A: Cross-candidate comparison of key stability metrics.",
            s_caption,
        ))
        story.append(Spacer(1, 8))

    cohort_data = {
        "Mean Stability": sum(a.stability.therapeutic_stability_score for a in analyses) / len(analyses),
        "Mean GC%": sum(a.stability.gc_content for a in analyses) / len(analyses) * 100,
        "Mean Confidence": sum(a.structure.structural_confidence for a in analyses) / len(analyses) * 100,
        "Mean Pairing": sum(a.structure.paired_fraction for a in analyses) / len(analyses) * 100,
    }
    story.append(CompBarChart(
        cohort_data, title="Cohort Average Metrics (%)", width=460, height=140,
    ))
    story.append(Spacer(1, 4))
    story.append(Paragraph("Figure B: Mean thermodynamic and structural metrics across the candidate cohort.", s_caption))
    story.append(Spacer(1, 8))

    story.append(PageBreak())

    # ================================================================== #
    # SECTION 3 — CANDIDATE DOSSIERS
    # ================================================================== #
    for idx, a in enumerate(analyses):
        stem_metrics = compute_stem_metrics(a.structure.dot_bracket, a.structure.base_pairs)

        story.append(Paragraph(f"3.{idx + 1} Candidate Dossier: {a.sequence_name}", s_section))
        story.append(HRFlowable(width="100%", thickness=0.5, color=PURPLE, spaceAfter=6))

        # -- 3.x.1 Therapeutic Classification ------------------------- #
        story.append(Paragraph(f"3.{idx + 1}.1 Therapeutic Classification &amp; Identification", s_subsection))

        clf_data = [
            [Paragraph("<b>Analysis ID</b>", s_table_body),
             Paragraph(a.analysis_id, s_table_body),
             Paragraph("<b>Classification</b>", s_table_body),
             Paragraph(
                 "High Viability" if a.stability.therapeutic_stability_score >= 70
                 else "Moderate Viability" if a.stability.therapeutic_stability_score >= 45
                 else "Limited Viability",
                 s_table_body)],
            [Paragraph("<b>Sequence Length</b>", s_table_body),
             Paragraph(f"{a.length} nt", s_table_body),
             Paragraph("<b>Molecular Type</b>", s_table_body),
             Paragraph("Single-Stranded RNA" if a.length > 50 else "Oligonucleotide", s_table_body)],
            [Paragraph("<b>GC Content</b>", s_table_body),
             Paragraph(f"{a.stability.gc_content * 100:.1f}%", s_table_body),
             Paragraph("<b>Stability Tier</b>", s_table_body),
             Paragraph(
                 "Tier I" if a.stability.therapeutic_stability_score >= 70
                 else "Tier II" if a.stability.therapeutic_stability_score >= 45
                 else "Tier III",
                 s_table_body)],
        ]
        t_clf = Table(clf_data, colWidths=[90, 140, 90, 140])
        t_clf.setStyle(TableStyle([
            ("ALIGN", (0, 0), (-1, -1), "LEFT"),
            ("TOPPADDING", (0, 0), (-1, -1), 4),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ("LINEBELOW", (0, 0), (-1, -1), 0.3, BORDER_GRAY),
        ]))
        story.append(t_clf)
        story.append(Spacer(1, 4))

        # Gauge meters
        gauge_data = [[
            GaugeMeter("Therapeutic Score", a.stability.therapeutic_stability_score, "/100", width=150, height=70, color=PURPLE),
            GaugeMeter("Structural Confidence", a.structure.structural_confidence * 100, "%", width=150, height=70, color=ACCENT_CYAN),
            GaugeMeter("GC Content", a.stability.gc_content * 100, "%", width=150, height=70, color=ACCENT_EMERALD),
        ]]
        t_gauges = Table(gauge_data, colWidths=[153, 153, 154])
        t_gauges.setStyle(TableStyle([
            ("ALIGN", (0, 0), (-1, -1), "CENTER"),
            ("LEFTPADDING", (0, 0), (-1, -1), 2),
            ("RIGHTPADDING", (0, 0), (-1, -1), 2),
        ]))
        story.append(t_gauges)
        story.append(Spacer(1, 6))

        # -- 3.x.2 Thermodynamic Analysis ----------------------------- #
        story.append(Paragraph(f"3.{idx + 1}.2 Thermodynamic Analysis", s_subsection))

        story.append(_build_metrics_table(a))
        story.append(Spacer(1, 6))

        story.append(GaugeBar("Therapeutic Stability Index", a.stability.therapeutic_stability_score, 100, width=460, height=20, color=PURPLE))
        story.append(Spacer(1, 2))
        story.append(GaugeBar("Structural Confidence", a.structure.structural_confidence * 100, 100, width=460, height=20, color=ACCENT_CYAN))
        story.append(Spacer(1, 2))
        story.append(GaugeBar("GC Content", a.stability.gc_content * 100, 100, width=460, height=20, color=ACCENT_EMERALD))
        story.append(Spacer(1, 2))
        story.append(GaugeBar("Paired Fraction", a.structure.paired_fraction * 100, 100, width=460, height=20, color=ACCENT_AMBER))
        story.append(Spacer(1, 8))

        thermo_text = _generate_thermo_interpretation(a)
        story.append(Paragraph(thermo_text, s_body))
        story.append(Spacer(1, 6))

        # Stem-loop metrics sub-section
        story.append(Paragraph(f"3.{idx + 1}.3 Stem-Loop Architecture Metrics", s_subsection))

        radar = RadarChart(
            categories=["Stability", "GC%", "Confidence", "Pairing", "Length\nScore"],
            values=[
                a.stability.therapeutic_stability_score,
                a.stability.gc_content * 100,
                a.structure.structural_confidence * 100,
                a.structure.paired_fraction * 100,
                min(100, a.length / 2),
            ],
            title="Stability Radar Profile",
            width=200, height=200,
        )

        comp_chart = CompBarChart(
            {k: v for k, v in [("A", a.composition.A * 100), ("C", a.composition.C * 100),
                               ("G", a.composition.G * 100), ("U", a.composition.U * 100)]},
            title="Nucleotide Distribution (%)", width=240, height=160,
        )

        radar_comp_data = [[radar, comp_chart]]
        t_radar = Table(radar_comp_data, colWidths=[210, 250])
        t_radar.setStyle(TableStyle([
            ("ALIGN", (0, 0), (-1, -1), "CENTER"),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ]))
        story.append(t_radar)
        story.append(Spacer(1, 4))

        stem_metrics_report = (
            f"<b>Stem Regions:</b> {stem_metrics['stem_count']} | "
            f"<b>Total Stem Pairs:</b> {stem_metrics['total_stem_pairs']} | "
            f"<b>Avg Stem Length:</b> {stem_metrics['avg_stem_length']} bp | "
            f"<b>Loop Regions:</b> {stem_metrics['loop_count']} | "
            f"<b>Avg Loop Size:</b> {stem_metrics['avg_loop_size']} nt | "
            f"<b>Max Pairing Depth:</b> {stem_metrics['max_depth']}"
        )
        story.append(Paragraph(stem_metrics_report, s_body_small))
        story.append(Spacer(1, 4))

        # Composition table
        story.append(Paragraph(f"3.{idx + 1}.4 Nucleotide Composition Analysis", s_subsection))
        story.append(_build_composition_table(a))
        story.append(Spacer(1, 6))

        # Sequence and dot-bracket
        story.append(Paragraph(f"<b>RNA Sequence (5' → 3'):</b>", s_body_small))
        seq_display = a.sequence
        if len(seq_display) > 100:
            seq_display = seq_display[:100] + "..."
        story.append(Paragraph(seq_display, s_code))

        story.append(Paragraph(f"<b>Dot-Bracket Secondary Fold Notation:</b>", s_body_small))
        db_display = a.structure.dot_bracket
        if len(db_display) > 130:
            db_display = db_display[:130] + "..."
        story.append(Paragraph(db_display, s_code))
        story.append(Spacer(1, 6))

        # -- 3.x.5 Structure Visualization ---------------------------- #
        story.append(Paragraph(f"3.{idx + 1}.5 RNA Secondary Structure Topology", s_subsection))
        story.append(Paragraph(
            "The following diagram presents the predicted secondary structure "
            "layout computed from dot-bracket notation, showing authentic stem-loop "
            "topology with color-coded nucleotide residues and base-pair connectivity.",
            s_body_small,
        ))
        story.append(Spacer(1, 4))

        fold_plot = StructureTopologyFlowable(
            sequence=a.sequence,
            dot_bracket=a.structure.dot_bracket,
            base_pairs=a.structure.base_pairs,
            width=460,
            height=260,
        )
        story.append(fold_plot)
        story.append(Spacer(1, 3))
        story.append(Paragraph(
            "Figure C: Predicted RNA secondary structure topology. Nucleotides colored by base type: "
            "A (green), C (blue), G (orange), U (red). Curved arcs represent base-pair interactions.",
            s_caption,
        ))
        story.append(Spacer(1, 6))

        # -- 3.x.6 Scientific Interpretation -------------------------- #
        story.append(Paragraph(f"3.{idx + 1}.6 Scientific Interpretation", s_subsection))

        interp_text = _generate_scientific_interpretation(a)
        sections = interp_text.split("\n\n")
        for sec in sections:
            if sec.strip():
                if sec.isupper():
                    story.append(Paragraph(f"<b>{sec}</b>", s_subsubsection))
                else:
                    story.append(Paragraph(sec, s_body))
        story.append(Spacer(1, 6))

        # -- 3.x.7 Risk Assessment ------------------------------------ #
        story.append(Paragraph(f"3.{idx + 1}.7 Risk Assessment", s_subsection))
        risks = _generate_risk_assessment(a)
        for risk_type, level, desc in risks:
            story.append(RiskBadge(risk_type, level, desc, width=460))
            story.append(Spacer(1, 2))

        story.append(Spacer(1, 8))

        # -- 3.x.8 Stability Verdict ---------------------------------- #
        story.append(Paragraph(f"3.{idx + 1}.8 Stability Verdict", s_subsection))
        score = a.stability.therapeutic_stability_score
        if score >= 75:
            verdict_text = (
                f"<b>High Viability</b> — Candidate '{a.sequence_name}' exhibits "
                f"robust thermodynamic properties and structural determinacy "
                f"suitable for therapeutic development."
            )
        elif score >= 50:
            verdict_text = (
                f"<b>Moderate Viability</b> — Candidate '{a.sequence_name}' demonstrates "
                f"acceptable stability characteristics; sequence optimization "
                f"may further enhance structural persistence."
            )
        elif score >= 25:
            verdict_text = (
                f"<b>Limited Viability</b> — Candidate '{a.sequence_name}' shows "
                f"marginal structural stability; engineering modifications "
                f"recommended to improve folding thermodynamics."
            )
        else:
            verdict_text = (
                f"<b>Low Viability</b> — Candidate '{a.sequence_name}' exhibits "
                f"insufficient structural stability for direct therapeutic "
                f"application in its current form."
            )
        story.append(Paragraph(verdict_text, s_body))
        story.append(Paragraph(
            f"Composite score: {score:.2f}/100 | "
            f"GC content: {a.stability.gc_content * 100:.1f}% | "
            f"MFE: {a.structure.minimum_free_energy:.2f} kcal/mol | "
            f"MFE/nt: {a.stability.mfe_per_nucleotide:.4f} kcal/mol/nt | "
            f"Ensemble diversity: {a.stability.ensemble_diversity:.4f}",
            s_body_small,
        ))
        story.append(Spacer(1, 12))

        if idx < len(analyses) - 1:
            story.append(PageBreak())

    # ================================================================== #
    # SECTION 4 — COMPARATIVE ANALYTICS
    # ================================================================== #
    if len(analyses) >= 2:
        story.append(PageBreak())
        story.append(Paragraph("4. Comparative Analytics", s_section))
        story.append(HRFlowable(width="100%", thickness=0.5, color=PURPLE, spaceAfter=8))

        sorted_analyses = sorted(
            analyses, key=lambda a: a.stability.therapeutic_stability_score, reverse=True
        )

        comp_headers = [
            Paragraph("<b>Rank</b>", s_table_header),
            Paragraph("<b>Candidate</b>", s_table_header),
            Paragraph("<b>Length</b>", s_table_header),
            Paragraph("<b>GC%</b>", s_table_header),
            Paragraph("<b>MFE</b>", s_table_header),
            Paragraph("<b>Paired%</b>", s_table_header),
            Paragraph("<b>Stems</b>", s_table_header),
            Paragraph("<b>Stability</b>", s_table_header),
        ]
        comp_rows = [comp_headers]
        for rank_idx, ca in enumerate(sorted_analyses):
            rank_label = f"#{rank_idx + 1}"
            if rank_idx == 0:
                rank_label = "★ #1"
            s_metrics = compute_stem_metrics(ca.structure.dot_bracket, ca.structure.base_pairs)
            comp_rows.append([
                Paragraph(rank_label, s_table_body),
                Paragraph(ca.sequence_name, s_table_body),
                Paragraph(f"{ca.length} nt", s_table_body),
                Paragraph(f"{ca.stability.gc_content * 100:.1f}%", s_table_body),
                Paragraph(f"{ca.structure.minimum_free_energy:.2f}", s_table_body),
                Paragraph(f"{ca.structure.paired_fraction * 100:.1f}%", s_table_body),
                Paragraph(str(s_metrics["stem_count"]), s_table_body),
                Paragraph(f"<b>{ca.stability.therapeutic_stability_score:.2f}</b>", s_table_body),
            ])

        num_cols = len(comp_headers)
        t_comp = Table(comp_rows, colWidths=[32, 100, 42, 38, 50, 48, 36, 54])
        t_comp.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), PURPLE),
            ("ALIGN", (0, 0), (-1, -1), "CENTER"),
            ("ALIGN", (1, 0), (1, -1), "LEFT"),
            ("TOPPADDING", (0, 0), (-1, -1), 5),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
            ("LINEBELOW", (0, 0), (-1, -1), 0.3, BORDER_GRAY),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1),
             [colors.HexColor("#0D1117"), colors.HexColor("#0A0E17")]),
        ]))
        story.append(t_comp)
        story.append(Spacer(1, 8))

        top_candidate = sorted_analyses[0]
        bottom_candidate = sorted_analyses[-1]
        spread = (
            top_candidate.stability.therapeutic_stability_score
            - bottom_candidate.stability.therapeutic_stability_score
        )
        story.append(Paragraph(
            f"<b>Cohort Analysis:</b> The {len(analyses)}-candidate cohort "
            f"spans a Therapeutic Stability Index range of "
            f"{bottom_candidate.stability.therapeutic_stability_score:.1f} "
            f"to {top_candidate.stability.therapeutic_stability_score:.1f} "
            f"(spread: {spread:.1f} points). "
            f"Candidate '{top_candidate.sequence_name}' leads the ranking "
            f"with a score of {top_candidate.stability.therapeutic_stability_score:.1f}/100.",
            s_body,
        ))
        story.append(Spacer(1, 6))

        story.append(Paragraph("Structural Relationship Summary", s_subsection))
        story.append(Paragraph(
            f"The comparative analysis encompasses {len(analyses)} candidate "
            f"sequences with diverse thermodynamic and structural profiles. "
            f"Molecular similarity assessment and structural clustering are "
            f"available through the platform's Molecular Similarity Desk "
            f"for detailed relationship mapping and family detection.",
            s_body_small,
        ))
        story.append(Spacer(1, 8))

    # ================================================================== #
    # SECTION 5 — TOP CANDIDATE SPOTLIGHT
    # ================================================================== #
    if len(analyses) >= 1:
        story.append(PageBreak())
        story.append(Paragraph("5. Top Candidate Spotlight", s_section))
        story.append(HRFlowable(width="100%", thickness=0.5, color=PURPLE, spaceAfter=8))

        top = max(analyses, key=lambda a: _safe_val(a.stability.therapeutic_stability_score))

        story.append(Paragraph(
            f"Lead Candidate: {top.sequence_name}",
            s_subsection,
        ))

        spotlight_data = [
            [Paragraph("<b>Metric</b>", s_table_header), Paragraph("<b>Value</b>", s_table_header),
             Paragraph("<b>Metric</b>", s_table_header), Paragraph("<b>Value</b>", s_table_header)],
            [Paragraph("Therapeutic Stability Index", s_table_body),
             Paragraph(f"{top.stability.therapeutic_stability_score:.2f}/100", s_table_body),
             Paragraph("GC Content", s_table_body),
             Paragraph(f"{top.stability.gc_content*100:.1f}%", s_table_body)],
            [Paragraph("Minimum Free Energy", s_table_body),
             Paragraph(f"{top.structure.minimum_free_energy:.2f} kcal/mol", s_table_body),
             Paragraph("MFE per Nucleotide", s_table_body),
             Paragraph(f"{top.stability.mfe_per_nucleotide:.4f} kcal/mol/nt", s_table_body)],
            [Paragraph("Structural Confidence", s_table_body),
             Paragraph(f"{top.structure.structural_confidence*100:.1f}%", s_table_body),
             Paragraph("Paired Fraction", s_table_body),
             Paragraph(f"{top.structure.paired_fraction*100:.1f}%", s_table_body)],
            [Paragraph("Ensemble Diversity", s_table_body),
             Paragraph(f"{top.stability.ensemble_diversity:.4f}", s_table_body),
             Paragraph("Sequence Length", s_table_body),
             Paragraph(f"{top.length} nt", s_table_body)],
        ]

        s_metrics = compute_stem_metrics(top.structure.dot_bracket, top.structure.base_pairs)
        spotlight_data.append([
            Paragraph("Stem Regions", s_table_body),
            Paragraph(str(s_metrics["stem_count"]), s_table_body),
            Paragraph("Loop Regions", s_table_body),
            Paragraph(str(s_metrics["loop_count"]), s_table_body),
        ])

        t_spot = Table(spotlight_data, colWidths=[110, 120, 110, 120])
        t_spot.setStyle(TableStyle([
            ("ALIGN", (0, 0), (-1, -1), "LEFT"),
            ("TOPPADDING", (0, 0), (-1, -1), 4),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ("GRID", (0, 0), (-1, -1), 0.3, BORDER_GRAY),
            ("BACKGROUND", (0, 0), (-1, 0), PURPLE),
            ("BACKGROUND", (1, 1), (1, -1), colors.HexColor("#0A0E17")),
            ("BACKGROUND", (3, 1), (3, -1), colors.HexColor("#0A0E17")),
            ("BACKGROUND", (0, 1), (0, -1), colors.HexColor("#0D1117")),
            ("BACKGROUND", (2, 1), (2, -1), colors.HexColor("#0D1117")),
        ]))
        story.append(t_spot)
        story.append(Spacer(1, 8))

        # Gauge bars for top candidate
        story.append(GaugeBar("Therapeutic Stability Index", top.stability.therapeutic_stability_score, 100, width=460, height=20, color=PURPLE))
        story.append(Spacer(1, 2))
        story.append(GaugeBar("GC Content", top.stability.gc_content * 100, 100, width=460, height=20, color=ACCENT_CYAN))
        story.append(Spacer(1, 2))
        story.append(GaugeBar("Structural Confidence", top.structure.structural_confidence * 100, 100, width=460, height=20, color=ACCENT_EMERALD))
        story.append(Spacer(1, 2))
        story.append(GaugeBar("Base Pairing Density", top.structure.paired_fraction * 100, 100, width=460, height=20, color=ACCENT_AMBER))
        story.append(Spacer(1, 8))

        # Scientific interpretation for top candidate
        try:
            interp = generate_scientific_interpretation(top)
            interp_paras = interp.split("\n\n")
            for para in interp_paras:
                if para.strip():
                    story.append(Paragraph(para, s_body))
                    story.append(Spacer(1, 4))
        except Exception:
            # Fallback to existing interpretation
            story.append(Paragraph(top.genomic_interpretation, s_body))

    # ================================================================== #
    # SECTION 6 — FINAL SCIENTIFIC CONCLUSION
    # ================================================================== #
    story.append(PageBreak())
    story.append(Paragraph("6. Scientific Conclusion", s_section))
    story.append(HRFlowable(width="100%", thickness=0.5, color=PURPLE, spaceAfter=8))

    if len(analyses) >= 2:
        try:
            comp = generate_comparative_summary(analyses)
            story.append(Paragraph(
                f"This dossier has presented a comprehensive computational assessment of "
                f"<b>{len(analyses)} RNA candidate(s)</b> for therapeutic development. "
                f"The cohort mean Therapeutic Stability Index of <b>{comp['mean_score']:.1f}/100</b> "
                f"reflects the integrated thermodynamic, structural, and compositional profile "
                f"of the candidate library.",
                s_body,
            ))
            story.append(Spacer(1, 4))

            story.append(Paragraph("Summary of Findings", s_subsubsection))
            story.append(Paragraph(comp["summary"], s_body))
            story.append(Spacer(1, 4))

            if comp["strengths"]:
                story.append(Paragraph("Cohort Strengths", s_subsubsection))
                for s in comp["strengths"]:
                    story.append(Paragraph(f"• {s}", s_body_small))
                story.append(Spacer(1, 4))

            if comp["weaknesses"]:
                story.append(Paragraph("Cohort Weaknesses", s_subsubsection))
                for w in comp["weaknesses"]:
                    story.append(Paragraph(f"• {w}", s_body_small))
                story.append(Spacer(1, 4))

            story.append(Paragraph("Recommendation", s_subsubsection))
            story.append(Paragraph(comp["recommendation"], s_body))
        except Exception:
            story.append(Paragraph(
                f"A total of <b>{len(analyses)} candidate(s)</b> were evaluated using "
                f"integrated thermodynamic, structural, and compositional analysis. "
                f"The top-ranked candidate demonstrates the most favorable balance of "
                f"structural stability, GC content, and conformational determinacy.",
                s_body,
            ))
    else:
        a = analyses[0]
        try:
            interp = generate_scientific_interpretation(a)
            paras = interp.split("\n\n")
            story.append(Paragraph(
                f"This dossier presents the computational assessment of "
                f"candidate <b>'{a.sequence_name}'</b>. The following conclusions "
                f"are drawn from the integrated thermodynamic and structural analysis.",
                s_body,
            ))
            story.append(Spacer(1, 4))
            for para in paras[:3]:
                if para.strip():
                    # Extract GC assessment, thermodynamic, compactness
                    story.append(Paragraph(para, s_body_small))
                    story.append(Spacer(1, 3))
        except Exception:
            story.append(Paragraph(
                f"Candidate '{a.sequence_name}' has been evaluated using the "
                f"RNA Genomics Intelligence Platform computational pipeline. "
                f"The analysis encompasses secondary structure prediction, "
                f"thermodynamic stability scoring, and therapeutic viability assessment.",
                s_body,
            ))

    story.append(Spacer(1, 10))
    story.append(HRFlowable(width="100%", thickness=0.5, color=BORDER_GRAY, spaceAfter=6))
    story.append(Paragraph(
        f"This concludes the computational genomics dossier. "
        f"All analyses were performed using the RNA Genomics Intelligence Platform. "
        f"The findings are intended to guide therapeutic candidate selection and "
        f"sequence optimization efforts.",
        s_caption,
    ))

    # ================================================================== #
    # SECTION 7 — EXPORT APPENDIX
    # ================================================================== #
    story.append(PageBreak())
    story.append(Paragraph("7. Export Appendix", s_section))
    story.append(HRFlowable(width="100%", thickness=0.5, color=PURPLE, spaceAfter=8))
    story.append(Paragraph(
        "The following appendix provides machine-readable export data "
        "and metadata for all candidates included in this dossier.",
        s_body_small,
    ))
    story.append(Spacer(1, 6))

    story.append(Paragraph("7.1 FASTA Nucleotide Sequences", s_subsection))
    for a in analyses:
        desc = (
            f"length={a.length} GC={a.stability.gc_content * 100:.1f}% "
            f"stability={a.stability.therapeutic_stability_score:.2f}/100 "
            f"MFE={a.structure.minimum_free_energy:.2f}"
        )
        fasta_entry = f">{a.sequence_name} {desc}\n{a.sequence}"
        story.append(Paragraph(fasta_entry, s_code))
        story.append(Spacer(1, 3))
    story.append(Spacer(1, 8))

    story.append(Paragraph("7.2 Analysis Metadata", s_subsection))
    meta_headers = [
        Paragraph("<b>Analysis ID</b>", s_table_header),
        Paragraph("<b>Candidate</b>", s_table_header),
        Paragraph("<b>Timestamp</b>", s_table_header),
        Paragraph("<b>Stability</b>", s_table_header),
        Paragraph("<b>MFE</b>", s_table_header),
    ]
    meta_rows = [meta_headers]
    for a in analyses:
        ts = a.timestamp[:19] if len(a.timestamp) > 19 else a.timestamp
        meta_rows.append([
            Paragraph(a.analysis_id, s_table_body),
            Paragraph(a.sequence_name, s_table_body),
            Paragraph(ts, s_table_body),
            Paragraph(f"{a.stability.therapeutic_stability_score:.2f}", s_table_body),
            Paragraph(f"{a.structure.minimum_free_energy:.2f}", s_table_body),
        ])

    t_meta2 = Table(meta_rows, colWidths=[90, 100, 120, 70, 60])
    t_meta2.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), PURPLE),
        ("ALIGN", (0, 0), (-1, -1), "LEFT"),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LINEBELOW", (0, 0), (-1, -1), 0.3, BORDER_GRAY),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1),
         [colors.HexColor("#0D1117"), colors.HexColor("#0A0E17")]),
    ]))
    story.append(t_meta2)
    story.append(Spacer(1, 8))

    story.append(Paragraph("7.3 Secondary Structure Notations", s_subsection))
    for a in analyses:
        story.append(Paragraph(
            f"<b>{a.sequence_name}</b> ({a.length} nt, "
            f"{len(a.structure.base_pairs)} bp, "
            f"{compute_stem_metrics(a.structure.dot_bracket, a.structure.base_pairs)['stem_count']} stems)",
            s_body_small,
        ))
        db = a.structure.dot_bracket
        if len(db) > 130:
            db = db[:130] + "..."
        story.append(Paragraph(db, s_code))
    story.append(Spacer(1, 10))

    story.append(HRFlowable(width="100%", thickness=0.5, color=BORDER_GRAY, spaceAfter=6))
    story.append(Paragraph(
        f"Generated by RNA Genomics Intelligence Platform | "
        f"{datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')} | "
        f"{len(analyses)} candidate(s) | END OF DOSSIER",
        s_caption,
    ))

    doc.build(story)
    buffer.seek(0)
    return buffer.getvalue()
