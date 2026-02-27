#!/usr/bin/env python3
"""
SLT Framework Architecture Diagram (Figure 1)

Creates a publication-ready schematic of the SLT classification cascade:
  - 4 evidence layers (POPAF, gnomAD, GERMQ, COSMIC)
  - Evidence level assignment (high/medium/low)
  - PureCN posterior thresholds (0.8/0.5/0.2)
  - CHIP classification module
  - 4-tier cascade (SLT-A → SLT-D)

Usage:
    python3 fig_slt_architecture.py --output figures/fig_slt_architecture.pdf
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np


# Color scheme
TIER_COLORS = {
    "SLT-A": "#1b9e77",
    "SLT-B": "#d95f02",
    "SLT-C": "#7570b3",
    "SLT-D": "#e7298a",
}

LAYER_COLORS = ["#a6cee3", "#b2df8a", "#fdbf6f", "#fb9a99"]
INPUT_COLOR = "#e0e0e0"
POSTERIOR_COLOR = "#cab2d6"
CHIP_COLOR = "#ffff99"
EVIDENCE_COLOR = "#f0f0f0"


def draw_rounded_box(ax, x, y, w, h, text, color, fontsize=8, fontweight="normal",
                     edgecolor="black", linewidth=1, text_color="black", alpha=1.0):
    """Draw a rounded rectangle with centered text."""
    box = FancyBboxPatch(
        (x - w/2, y - h/2), w, h,
        boxstyle="round,pad=0.02",
        facecolor=color, edgecolor=edgecolor, linewidth=linewidth,
        alpha=alpha,
    )
    ax.add_patch(box)
    ax.text(x, y, text, ha="center", va="center",
            fontsize=fontsize, fontweight=fontweight, color=text_color)


def draw_arrow(ax, x1, y1, x2, y2, color="black", style="-|>", lw=1.2):
    """Draw an arrow between two points."""
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle=style, color=color, lw=lw))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    ax.set_xlim(-0.5, 11.5)
    ax.set_ylim(-0.5, 8.5)
    ax.axis("off")

    # Title
    ax.text(5.5, 8.2, "Somatic Likelihood Tiering (SLT) — Classification Framework",
            ha="center", fontsize=14, fontweight="bold")

    # ─── INPUT BOX ───────────────────────────────────────────────────────────
    draw_rounded_box(ax, 2.5, 7.3, 4.0, 0.6,
                     "Tumor-Only WES Variants\n(Mutect2 unfiltered VCF)",
                     INPUT_COLOR, fontsize=9, fontweight="bold")

    draw_arrow(ax, 2.5, 7.0, 2.5, 6.6)

    # ─── EVIDENCE LAYERS (LEFT COLUMN) ───────────────────────────────────────
    ax.text(2.5, 6.5, "Evidence Layers", ha="center", fontsize=10,
            fontweight="bold", style="italic")

    layers = [
        ("Layer 1: POPAF", "POPAF ≥ 5.0\n(somatic-supporting)", LAYER_COLORS[0]),
        ("Layer 2: gnomAD", "gnomAD AF < 0.001\nor missing", LAYER_COLORS[1]),
        ("Layer 3: GERMQ", "GERMQ ≥ 30\n(germline quality)", LAYER_COLORS[2]),
        ("Layer 4: COSMIC", "Confirmed ≥ 5 OR\nhotspot ≥ 10 OR CGC", LAYER_COLORS[3]),
    ]

    y_start = 5.9
    for i, (title, desc, color) in enumerate(layers):
        y = y_start - i * 0.85
        draw_rounded_box(ax, 2.5, y, 3.8, 0.7,
                         f"{title}\n{desc}", color, fontsize=7)

    # Arrow from layers to evidence level
    draw_arrow(ax, 4.5, 4.2, 5.8, 4.2, style="-|>")

    # ─── EVIDENCE LEVEL ──────────────────────────────────────────────────────
    draw_rounded_box(ax, 7.2, 4.2, 2.6, 1.2,
                     "Evidence Level\n\nhigh: ≥3 layers + COSMIC\nmedium: ≥2 layers\nlow: 0-1 layers",
                     EVIDENCE_COLOR, fontsize=7)

    # ─── PURECN POSTERIOR (RIGHT COLUMN) ──────────────────────────────────────
    draw_rounded_box(ax, 8.5, 7.3, 3.5, 0.6,
                     "PureCN Tumor-Only\n(POSTERIOR.SOMATIC probability)",
                     POSTERIOR_COLOR, fontsize=9, fontweight="bold")

    draw_arrow(ax, 8.5, 7.0, 8.5, 6.6)

    # Posterior thresholds
    ax.text(8.5, 6.5, "Posterior Thresholds", ha="center", fontsize=10,
            fontweight="bold", style="italic")

    thresholds = [
        ("≥ 0.80", "SLT-A gate", "#c7e9c0"),
        ("≥ 0.50", "SLT-B gate", "#fdd0a2"),
        ("≥ 0.20", "SLT-C gate", "#c6dbef"),
    ]

    for i, (thresh, desc, color) in enumerate(thresholds):
        y = 5.9 - i * 0.75
        draw_rounded_box(ax, 8.5, y, 2.8, 0.55,
                         f"POSTERIOR ≥ {thresh.split()[1]}  →  {desc}", color, fontsize=7.5)

    # ─── CHIP MODULE ─────────────────────────────────────────────────────────
    draw_rounded_box(ax, 2.5, 2.5, 3.8, 1.0,
                     "CHIP Classification\n\nTier 1: 13 genes (DNMT3A, TET2, ...)\n"
                     "Tier 2: 28 genes (TP53, KRAS, ...)\n"
                     "Hotspot + VAF thresholds",
                     CHIP_COLOR, fontsize=7)

    draw_arrow(ax, 4.5, 2.5, 5.5, 1.5, style="-|>")

    # ─── CASCADE (BOTTOM) ────────────────────────────────────────────────────
    ax.text(5.5, 1.9, "SLT Tier Assignment Cascade", ha="center", fontsize=10,
            fontweight="bold", style="italic")

    draw_arrow(ax, 7.2, 3.5, 7.2, 1.5, style="-|>")
    draw_arrow(ax, 8.5, 4.3, 8.5, 1.5, style="-|>")

    # Tier boxes
    tier_info = [
        ("SLT-A", "High confidence\nP ≥ 0.8 + evidence\n+ not CHIP"),
        ("SLT-B", "Likely somatic\nP ≥ 0.5 OR\nhigh evidence"),
        ("SLT-C", "Possible somatic\nP ≥ 0.2 OR\nmedium evidence"),
        ("SLT-D", "Unlikely somatic\nEverything else"),
    ]

    x_start = 2.0
    for i, (tier, desc) in enumerate(tier_info):
        x = x_start + i * 2.5
        draw_rounded_box(ax, x, 0.6, 2.2, 1.2,
                         f"{tier}\n\n{desc}",
                         TIER_COLORS[tier], fontsize=7.5, fontweight="bold",
                         text_color="white", linewidth=1.5)

    # Arrow between tiers showing cascade
    for i in range(3):
        x = x_start + i * 2.5 + 1.2
        ax.annotate("", xy=(x + 0.2, 0.6), xytext=(x - 0.1, 0.6),
                    arrowprops=dict(arrowstyle="->", color="gray", lw=1))
        ax.text(x + 0.05, 0.1, "else", fontsize=6, ha="center", color="gray")

    # ─── Annotations ─────────────────────────────────────────────────────────
    ax.text(0.0, -0.3,
            "SLT cascade: first matching condition wins. "
            "Degraded mode (no PureCN): Layers 1,3 disabled, posterior = 0.",
            fontsize=7, color="gray", style="italic")

    # Save
    outpath = Path(args.output)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(outpath), dpi=300, bbox_inches="tight")
    fig.savefig(str(outpath.with_suffix(".png")), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")
    print(f"Saved: {outpath.with_suffix('.png')}")


if __name__ == "__main__":
    main()
