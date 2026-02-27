#!/usr/bin/env python3
"""
MC3 BLCA Degraded-Mode SLT Analysis Figure

Creates a 3-panel figure showing SLT degraded mode applied to TCGA MC3 BLCA:
  (A) Tier distribution pie/bar chart
  (B) Top genes in SLT-C tier
  (C) Variant classification breakdown for SLT-C

Usage:
    python3 fig_mc3_degraded.py \
        --slt-tsv work/mc3/mc3_slt_classified.tsv \
        --output  figures/fig_mc3_degraded.pdf
"""

import argparse
import csv
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


TIER_COLORS = {
    "SLT-A": "#1b9e77",
    "SLT-B": "#d95f02",
    "SLT-C": "#7570b3",
    "SLT-D": "#e7298a",
}

VARCLASS_COLORS = {
    "Missense_Mutation": "#4393c3",
    "Nonsense_Mutation": "#d6604d",
    "Frame_Shift_Del": "#f4a582",
    "Frame_Shift_Ins": "#fddbc7",
    "Splice_Site": "#92c5de",
    "Silent": "#969696",
    "In_Frame_Del": "#b2abd2",
    "In_Frame_Ins": "#e7d4e8",
    "Other": "#d9d9d9",
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--slt-tsv", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Parse classified TSV
    tiers = Counter()
    sltc_genes = Counter()
    sltc_varclass = Counter()
    samples = set()
    sltc_samples = set()

    with open(args.slt_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            tier = row["slt_tier"]
            tiers[tier] += 1
            sample = row["Tumor_Sample_Barcode"]
            samples.add(sample)

            if tier == "SLT-C":
                sltc_genes[row["Hugo_Symbol"]] += 1
                vc = row["Variant_Classification"]
                sltc_varclass[vc] += 1
                sltc_samples.add(sample)

    total = sum(tiers.values())
    n_samples = len(samples)
    n_sltc_samples = len(sltc_samples)

    print(f"Total variants: {total}")
    print(f"Total samples: {n_samples}")
    print(f"Tier distribution: {dict(tiers)}")
    print(f"SLT-C samples: {n_sltc_samples}/{n_samples}")

    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5),
                              gridspec_kw={"width_ratios": [1.0, 1.5, 1.0]})

    # ─── Panel A: Tier distribution ─────────────────────────────────────────
    ax = axes[0]
    tier_order = ["SLT-A", "SLT-B", "SLT-C", "SLT-D"]
    counts = [tiers.get(t, 0) for t in tier_order]
    pcts = [100 * c / total for c in counts]

    bars = ax.bar(tier_order, counts,
                  color=[TIER_COLORS[t] for t in tier_order],
                  edgecolor="black", linewidth=0.5)

    for bar, count, pct in zip(bars, counts, pcts):
        if count > 0:
            y_pos = bar.get_height()
            label = f"{count:,}\n({pct:.1f}%)" if pct >= 0.1 else f"{count:,}\n(<0.1%)"
            ax.text(bar.get_x() + bar.get_width()/2, y_pos,
                    label, ha="center", va="bottom", fontsize=7)

    ax.set_ylabel("Variant count", fontsize=9)
    ax.set_title("(A) SLT Tier Distribution\n(Degraded Mode)", fontsize=10, fontweight="bold")
    ax.set_yscale("log")
    ax.set_ylim(1, total * 5)
    ax.tick_params(labelsize=8)

    # Add summary text
    ax.text(0.98, 0.98,
            f"N = {total:,} variants\n{n_samples} BLCA samples\nDegraded mode\n(no PureCN)",
            transform=ax.transAxes, fontsize=7, va="top", ha="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    # ─── Panel B: Top genes in SLT-C ────────────────────────────────────────
    ax = axes[1]
    top_n = 15
    top_genes = sltc_genes.most_common(top_n)
    genes = [g for g, _ in top_genes]
    gene_counts = [c for _, c in top_genes]

    y_pos = np.arange(len(genes))
    bars = ax.barh(y_pos, gene_counts, color=TIER_COLORS["SLT-C"],
                   edgecolor="black", linewidth=0.5, alpha=0.85)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(genes, fontsize=8, fontstyle="italic")
    ax.invert_yaxis()
    ax.set_xlabel("Variant count (SLT-C)", fontsize=9)
    ax.set_title("(B) Top Genes in SLT-C Tier", fontsize=10, fontweight="bold")
    ax.tick_params(labelsize=8)

    # Add count labels
    for bar, count in zip(bars, gene_counts):
        ax.text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2,
                str(count), ha="left", va="center", fontsize=7)

    # Highlight known bladder cancer drivers
    known_drivers = {"TP53", "PIK3CA", "CDKN2A", "FGFR3", "KRAS", "HRAS",
                     "NFE2L2", "FBXW7", "PTEN", "CTNNB1", "ERBB2", "EGFR"}
    for i, gene in enumerate(genes):
        if gene in known_drivers:
            ax.get_yticklabels()[i].set_fontweight("bold")
            ax.get_yticklabels()[i].set_color("#1a1a1a")

    ax.text(0.98, 0.98,
            f"{len(sltc_genes)} genes total\n{n_sltc_samples}/{n_samples} samples\nwith SLT-C variants",
            transform=ax.transAxes, fontsize=7, va="top", ha="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    # ─── Panel C: Variant classification pie ─────────────────────────────────
    ax = axes[2]

    # Group minor categories
    main_classes = ["Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
                    "Frame_Shift_Ins", "Splice_Site", "Silent"]
    vc_data = {}
    other_count = 0
    for vc, count in sltc_varclass.items():
        if vc in main_classes:
            vc_data[vc] = count
        else:
            other_count += count
    if other_count > 0:
        vc_data["Other"] = other_count

    labels = list(vc_data.keys())
    sizes = list(vc_data.values())
    colors = [VARCLASS_COLORS.get(l, "#d9d9d9") for l in labels]

    # Clean up labels for display
    display_labels = [l.replace("_", " ").replace("Mutation", "").strip() for l in labels]

    wedges, texts, autotexts = ax.pie(
        sizes, labels=display_labels, colors=colors,
        autopct=lambda pct: f"{pct:.0f}%" if pct >= 3 else "",
        startangle=90, textprops={"fontsize": 7},
        wedgeprops={"edgecolor": "white", "linewidth": 0.5}
    )
    for autotext in autotexts:
        autotext.set_fontsize(7)

    ax.set_title("(C) Variant Classification\n(SLT-C Only)", fontsize=10, fontweight="bold")

    fig.tight_layout(w_pad=2.0)

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
