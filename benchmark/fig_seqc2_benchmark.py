#!/usr/bin/env python3
"""
SEQC2 HCC1395 Benchmark — SLT Performance Figures

Generates publication-ready multi-panel figure:
  A) Cumulative sensitivity–precision trade-off by SLT tier
  B) Per-tier variant classification (TP/FP composition)
  C) Variant type breakdown per tier (SNV vs indel)
  D) POSTERIOR.SOMATIC distribution for TP vs FP variants

Usage:
    python3 fig_seqc2_benchmark.py \
        --eval-tsv     slt/evaluation_per_variant.tsv \
        --slt-tsv      slt/slt_classified.tsv \
        --eval-bed     data/evaluation_regions/final_evaluation_regions_with_cov.bed \
        --eval-summary slt/evaluation_summary.txt \
        --output       figures/fig_seqc2_benchmark.pdf
"""

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np


# ── Style ────────────────────────────────────────────────────────────────────
TIER_COLORS = {
    "SLT-A": "#1b9e77",   # dark teal
    "SLT-B": "#d95f02",   # orange
    "SLT-C": "#7570b3",   # purple
    "SLT-D": "#e7298a",   # magenta
}
TIER_ORDER = ["SLT-A", "SLT-B", "SLT-C", "SLT-D"]
CUM_LABELS = ["≥SLT-A", "≥SLT-B", "≥SLT-C", "≥SLT-D"]

TP_COLOR = "#2166ac"
FP_COLOR = "#b2182b"
FN_COLOR = "#969696"


def parse_eval_tsv(path):
    """Parse per-variant evaluation TSV."""
    variants = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            variants.append(row)
    return variants


def parse_slt_classified(path):
    """Parse SLT-classified TSV for posterior distribution analysis."""
    variants = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            variants.append(row)
    return variants


def parse_eval_bed(path):
    """Parse evaluation BED to get regions set."""
    regions = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            regions.setdefault(chrom, []).append((start, end))
    # Sort regions
    for chrom in regions:
        regions[chrom].sort()
    return regions


def in_region(chrom, pos, regions):
    """Check if position falls in evaluation regions."""
    if chrom not in regions:
        return False
    import bisect
    intervals = regions[chrom]
    idx = bisect.bisect_right([iv[0] for iv in intervals], pos) - 1
    if idx >= 0 and intervals[idx][0] <= pos < intervals[idx][1]:
        return True
    # Also check idx+1 in case of ties
    if idx + 1 < len(intervals) and intervals[idx + 1][0] <= pos < intervals[idx + 1][1]:
        return True
    return False


def parse_eval_summary(path):
    """Parse evaluation_summary.txt to get authoritative truth count and metrics."""
    summary = {"n_truth": None, "metrics": {}}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("Truth in eval regions:"):
                summary["n_truth"] = int(line.split(":")[1].strip())
            elif line.startswith("Truth SNVs in eval regions:"):
                summary["n_truth_snv"] = int(line.split(":")[1].strip())
            # Parse cumulative metric rows from evaluate_seqc2.py output.
            for prefix in CUM_LABELS:
                if line.startswith(prefix) or line.startswith(prefix.replace("≥", ">=")):
                    parts = line.split()
                    if len(parts) >= 8:
                        summary["metrics"][prefix] = {
                            "called": int(parts[1]),
                            "tp": int(parts[2]),
                            "fp": int(parts[3]),
                            "fn": int(parts[4]),
                            "sens": float(parts[5]),
                            "prec": float(parts[6]),
                            "f1": float(parts[7]),
                        }
    return summary


def compute_cumulative_metrics(eval_variants, n_truth_override=None):
    """Compute sensitivity, precision, F1 at cumulative tier thresholds."""
    # Only use variants in eval regions
    in_eval = [v for v in eval_variants if v.get("in_eval_region", "True") == "True"]

    # Count truth variants that Mutect2 called (detectable TPs)
    truth_in_eval = set()
    for v in in_eval:
        if v["in_truth"] == "True":
            key = (v["chrom"], v["pos"], v["ref"], v["alt"])
            truth_in_eval.add(key)

    # Use authoritative truth count from evaluate_seqc2.py if provided
    n_truth = n_truth_override if n_truth_override else len(truth_in_eval)

    metrics = {}
    for i, tier in enumerate(TIER_ORDER):
        # Cumulative: include this tier and all higher tiers
        included_tiers = set(TIER_ORDER[:i + 1])

        called = [v for v in in_eval
                  if v["tier"] in included_tiers]
        tp = sum(1 for v in called if v["in_truth"] == "True")
        fp = sum(1 for v in called if v["in_truth"] == "False")
        fn = n_truth - tp

        sens = tp / n_truth if n_truth > 0 else 0
        prec = tp / (tp + fp) if (tp + fp) > 0 else 0
        f1 = 2 * sens * prec / (sens + prec) if (sens + prec) > 0 else 0

        metrics[tier] = {
            "called": len(called), "tp": tp, "fp": fp, "fn": fn,
            "sens": sens, "prec": prec, "f1": f1
        }

    return metrics, n_truth


def compute_per_tier_counts(eval_variants):
    """Count TP/FP per individual tier."""
    in_eval = [v for v in eval_variants if v.get("in_eval_region", "True") == "True"]
    counts = {}
    for tier in TIER_ORDER:
        tier_vars = [v for v in in_eval if v["tier"] == tier]
        tp = sum(1 for v in tier_vars if v["in_truth"] == "True")
        fp = sum(1 for v in tier_vars if v["in_truth"] == "False")
        counts[tier] = {"tp": tp, "fp": fp, "total": len(tier_vars)}
    return counts


def get_posterior_distributions(slt_variants, eval_variants, eval_regions):
    """Get POSTERIOR.SOMATIC distributions for TP vs FP variants in eval regions."""
    # Build truth set from eval variants
    truth_keys = set()
    for v in eval_variants:
        if v["in_truth"] == "True" and v.get("in_eval_region", "True") == "True":
            truth_keys.add((v["chrom"], v["pos"], v["ref"], v["alt"]))

    tp_posteriors = []
    fp_posteriors = []

    for v in slt_variants:
        chrom = v.get("Chromosome", "")
        pos = v.get("Start_Position", "")
        ref = v.get("Reference_Allele", "")
        alt = v.get("Tumor_Seq_Allele2", "")

        if not chrom or not pos:
            continue

        # Check if in eval region
        if not in_region(chrom, int(pos), eval_regions):
            continue

        posterior = v.get("POSTERIOR.SOMATIC", "")
        if not posterior or posterior == "NA" or posterior == "":
            continue
        try:
            posterior = float(posterior)
        except ValueError:
            continue

        key = (chrom, pos, ref, alt)
        if key in truth_keys:
            tp_posteriors.append(posterior)
        else:
            fp_posteriors.append(posterior)

    return tp_posteriors, fp_posteriors


def plot_figure(cum_metrics, per_tier, n_truth, tp_post, fp_post, outpath):
    """Create multi-panel benchmark figure."""
    fig, axes = plt.subplots(2, 2, figsize=(10, 8.5))
    fig.subplots_adjust(hspace=0.35, wspace=0.30, top=0.93, bottom=0.08,
                        left=0.10, right=0.95)

    # ── Panel A: Sensitivity–Precision trade-off ─────────────────────────────
    ax = axes[0, 0]
    sens_vals = [cum_metrics[t]["sens"] for t in TIER_ORDER]
    prec_vals = [cum_metrics[t]["prec"] for t in TIER_ORDER]
    f1_vals = [cum_metrics[t]["f1"] for t in TIER_ORDER]

    # Plot connected line
    ax.plot(sens_vals, prec_vals, "k-", linewidth=1.5, zorder=2)

    # Plot tier points
    for i, tier in enumerate(TIER_ORDER):
        ax.scatter(sens_vals[i], prec_vals[i], c=TIER_COLORS[tier],
                   s=120, zorder=3, edgecolors="black", linewidth=0.8)
        # Label each point
        offset_x = -0.05 if i >= 2 else 0.03
        offset_y = 0.03 if i < 2 else -0.05
        ax.annotate(CUM_LABELS[i],
                    (sens_vals[i], prec_vals[i]),
                    xytext=(sens_vals[i] + offset_x, prec_vals[i] + offset_y),
                    fontsize=8, fontweight="bold",
                    ha="center" if i >= 2 else "left")

    ax.set_xlabel("Sensitivity", fontsize=10)
    ax.set_ylabel("Precision", fontsize=10)
    ax.set_title("Cumulative Sensitivity–Precision", fontsize=11, fontweight="bold")
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.axhline(y=0.5, color="gray", linestyle=":", alpha=0.4)
    ax.axvline(x=0.5, color="gray", linestyle=":", alpha=0.4)

    # Add F1 isolines
    for f1_target in [0.1, 0.2, 0.3, 0.5]:
        s_range = np.linspace(0.01, 1.0, 200)
        p_range = (f1_target * s_range) / (2 * s_range - f1_target)
        valid = (p_range > 0) & (p_range <= 1)
        ax.plot(s_range[valid], p_range[valid], ":", color="#cccccc",
                linewidth=0.7, alpha=0.6)
        # Label the isoline
        idx = np.argmin(np.abs(s_range[valid] - 0.95))
        if valid.sum() > 0:
            ax.text(0.98, p_range[valid][min(idx, len(p_range[valid])-1)],
                    f"F1={f1_target}", fontsize=6, color="#999999",
                    ha="right", va="bottom")

    # ── Panel B: Per-tier TP/FP stacked bar ──────────────────────────────────
    ax = axes[0, 1]
    x = np.arange(len(TIER_ORDER))
    tp_vals = [per_tier[t]["tp"] for t in TIER_ORDER]
    fp_vals = [per_tier[t]["fp"] for t in TIER_ORDER]

    bars_tp = ax.bar(x, tp_vals, width=0.6, color=TP_COLOR, label="True Positive",
                     edgecolor="white", linewidth=0.5)
    bars_fp = ax.bar(x, fp_vals, width=0.6, bottom=tp_vals, color=FP_COLOR,
                     label="False Positive", edgecolor="white", linewidth=0.5)

    # Add count labels
    for i, tier in enumerate(TIER_ORDER):
        total = tp_vals[i] + fp_vals[i]
        if total > 0:
            prec = tp_vals[i] / total * 100
            ax.text(i, total + max(fp_vals) * 0.02,
                    f"{prec:.0f}%", ha="center", va="bottom", fontsize=8,
                    fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(TIER_ORDER, fontsize=9)
    ax.set_ylabel("Variant Count", fontsize=10)
    ax.set_title("Per-Tier Classification", fontsize=11, fontweight="bold")
    ax.legend(fontsize=8, loc="upper right")
    ax.set_ylim(0, max(tp_vals[i] + fp_vals[i] for i in range(len(TIER_ORDER))) * 1.15)

    # ── Panel C: Cumulative metrics bar chart ────────────────────────────────
    ax = axes[1, 0]
    x = np.arange(len(TIER_ORDER))
    width = 0.25

    sens_bars = [cum_metrics[t]["sens"] for t in TIER_ORDER]
    prec_bars = [cum_metrics[t]["prec"] for t in TIER_ORDER]
    f1_bars = [cum_metrics[t]["f1"] for t in TIER_ORDER]

    ax.bar(x - width, sens_bars, width, label="Sensitivity",
           color="#2166ac", edgecolor="white", linewidth=0.5)
    ax.bar(x, prec_bars, width, label="Precision",
           color="#b2182b", edgecolor="white", linewidth=0.5)
    ax.bar(x + width, f1_bars, width, label="F1 Score",
           color="#4daf4a", edgecolor="white", linewidth=0.5)

    # Add value labels
    for i in range(len(TIER_ORDER)):
        ax.text(i - width, sens_bars[i] + 0.02,
                f"{sens_bars[i]:.2f}", ha="center", va="bottom", fontsize=7)
        ax.text(i, prec_bars[i] + 0.02,
                f"{prec_bars[i]:.2f}", ha="center", va="bottom", fontsize=7)
        ax.text(i + width, f1_bars[i] + 0.02,
                f"{f1_bars[i]:.2f}", ha="center", va="bottom", fontsize=7)

    ax.set_xticks(x)
    ax.set_xticklabels(CUM_LABELS, fontsize=9)
    ax.set_ylabel("Score", fontsize=10)
    ax.set_title("Cumulative Tier Performance", fontsize=11, fontweight="bold")
    ax.legend(fontsize=8, loc="upper center", ncol=3)
    ax.set_ylim(0, 1.15)
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1.0))

    # ── Panel D: POSTERIOR.SOMATIC distribution (TP vs FP) ───────────────────
    ax = axes[1, 1]

    if tp_post and fp_post:
        bins = np.linspace(0, 1, 30)
        ax.hist(tp_post, bins=bins, alpha=0.7, color=TP_COLOR,
                label=f"TP (n={len(tp_post)})", density=True, edgecolor="white",
                linewidth=0.5)
        ax.hist(fp_post, bins=bins, alpha=0.7, color=FP_COLOR,
                label=f"FP (n={len(fp_post)})", density=True, edgecolor="white",
                linewidth=0.5)

        # Add tier threshold lines
        for thresh, label in [(0.8, "SLT-A"), (0.5, "SLT-B"), (0.2, "SLT-C")]:
            ax.axvline(x=thresh, color="gray", linestyle="--", linewidth=1, alpha=0.7)
            ax.text(thresh + 0.01, ax.get_ylim()[1] * 0.95, label,
                    fontsize=7, rotation=90, va="top", color="gray")

        ax.set_xlabel("POSTERIOR.SOMATIC", fontsize=10)
        ax.set_ylabel("Density", fontsize=10)
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "No PureCN posterior\ndata available",
                ha="center", va="center", fontsize=10, transform=ax.transAxes)

    ax.set_title("Posterior Distribution: TP vs FP", fontsize=11, fontweight="bold")

    # ── Panel labels ─────────────────────────────────────────────────────────
    for i, (ax, label) in enumerate(zip(axes.flat, ["A", "B", "C", "D"])):
        ax.text(-0.12, 1.08, label, transform=ax.transAxes,
                fontsize=14, fontweight="bold", va="top")

    # ── Annotation ───────────────────────────────────────────────────────────
    fig.text(0.5, 0.01,
             f"SEQC2 HCC1395 WES Benchmark | Truth: {n_truth} variants in evaluation regions | "
             f"Mutect2 tumor-only → SLT classification",
             ha="center", fontsize=8, color="gray")

    # Save
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(str(outpath), dpi=300, bbox_inches="tight")
    print(f"Saved: {outpath}")

    # Also save PNG
    png_path = outpath.with_suffix(".png")
    fig.savefig(str(png_path), dpi=300, bbox_inches="tight")
    print(f"Saved: {png_path}")

    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="SEQC2 SLT Benchmark Figures")
    parser.add_argument("--eval-tsv", required=True,
                        help="Per-variant evaluation TSV")
    parser.add_argument("--slt-tsv", required=True,
                        help="SLT-classified TSV")
    parser.add_argument("--eval-bed", required=True,
                        help="Evaluation regions BED")
    parser.add_argument("--eval-summary", required=True,
                        help="Evaluation summary TXT (from evaluate_seqc2.py)")
    parser.add_argument("--output", required=True,
                        help="Output figure path (PDF)")
    args = parser.parse_args()

    print("Loading evaluation data...")
    eval_variants = parse_eval_tsv(args.eval_tsv)
    print(f"  {len(eval_variants)} evaluated variants")

    print("Loading SLT classifications...")
    slt_variants = parse_slt_classified(args.slt_tsv)
    print(f"  {len(slt_variants)} SLT variants")

    print("Loading evaluation regions...")
    eval_regions = parse_eval_bed(args.eval_bed)
    n_regions = sum(len(v) for v in eval_regions.values())
    print(f"  {n_regions} regions")

    print("Loading evaluation summary...")
    summary = parse_eval_summary(args.eval_summary)
    n_truth_from_summary = summary.get("n_truth")
    print(f"  Truth in eval regions: {n_truth_from_summary}")

    print("Computing metrics...")
    cum_metrics, n_truth = compute_cumulative_metrics(
        eval_variants, n_truth_override=n_truth_from_summary)
    per_tier = compute_per_tier_counts(eval_variants)

    print("Extracting posterior distributions...")
    tp_post, fp_post = get_posterior_distributions(
        slt_variants, eval_variants, eval_regions)
    print(f"  TP with posterior: {len(tp_post)}, FP with posterior: {len(fp_post)}")

    print("Generating figure...")
    plot_figure(cum_metrics, per_tier, n_truth, tp_post, fp_post, args.output)

    # Print summary table
    print("\nSummary:")
    print(f"{'Threshold':<15} {'Called':>7} {'TP':>5} {'FP':>5} {'FN':>5} "
          f"{'Sens':>7} {'Prec':>7} {'F1':>7}")
    print("-" * 65)
    for i, tier in enumerate(TIER_ORDER):
        m = cum_metrics[tier]
        print(f"{CUM_LABELS[i]:<15} {m['called']:>7} {m['tp']:>5} {m['fp']:>5} "
              f"{m['fn']:>5} {m['sens']:>7.4f} {m['prec']:>7.4f} {m['f1']:>7.4f}")


if __name__ == "__main__":
    main()
