#!/usr/bin/env python3
"""SEQC2 Benchmark Evaluation — SLT vs Truth Set.

Protocol: 
Compares SLT-classified variants against SEQC2 HCC1395 truth set to compute
sensitivity, specificity, precision, and F1 at each SLT tier threshold.

Usage:
    python evaluate_seqc2.py \
        --slt-input slt_classified.tsv \
        --truth-snv truth_snv_eval_regions.vcf.gz \
        --truth-indel truth_indel_eval_regions.vcf.gz \
        --eval-bed final_evaluation_regions_with_cov.bed \
        --output evaluation_results.tsv \
        --summary evaluation_summary.txt
"""

import argparse
import csv
import gzip
import sys
from collections import defaultdict


def parse_vcf_truth(vcf_path):
    """Parse truth VCF into a set of (chrom, pos, ref, alt) tuples."""
    truth = set()
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            for alt in fields[4].split(","):
                truth.add((chrom, pos, ref, alt))
    return truth


def parse_bed_regions(bed_path):
    """Parse BED file into list of (chrom, start, end) tuples."""
    regions = []
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("track"):
                continue
            fields = line.strip().split("\t")
            if len(fields) >= 3:
                regions.append((fields[0], int(fields[1]), int(fields[2])))
    return regions


def variant_in_regions(chrom, pos, regions_by_chrom):
    """Check if a variant falls within evaluation regions (binary search)."""
    if chrom not in regions_by_chrom:
        return False
    for start, end in regions_by_chrom[chrom]:
        if start <= pos <= end:
            return True
    return False


def build_region_index(regions):
    """Index BED regions by chromosome for faster lookup."""
    by_chrom = defaultdict(list)
    for chrom, start, end in regions:
        by_chrom[chrom].append((start, end))
    # Sort by start for potential binary search
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    return by_chrom


def parse_slt_variants(slt_path):
    """Parse SLT classified TSV into list of variant dicts."""
    variants = []
    with open(slt_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            chrom = row.get("Chromosome") or row.get("chr") or row.get("CHROM", "")
            pos_str = row.get("Start_Position") or row.get("start") or row.get("POS", "")
            ref = row.get("Reference_Allele") or row.get("REF", "")
            alt = row.get("Tumor_Seq_Allele2") or row.get("ALT", "")
            tier = row.get("slt_tier", "")
            posterior = row.get("slt_posterior_used", "")
            gene = row.get("Hugo_Symbol", "")
            consequence = row.get("Variant_Classification", "")
            filt = row.get("FILTER", "")

            try:
                pos = int(pos_str)
            except (ValueError, TypeError):
                continue

            variants.append({
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "tier": tier,
                "posterior": posterior,
                "gene": gene,
                "consequence": consequence,
                "filter": filt,
                "key": (chrom, pos, ref, alt),
            })
    return variants


def compute_metrics(tp, fp, fn):
    """Compute sensitivity, precision, F1."""
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = 2 * precision * sensitivity / (precision + sensitivity) if (precision + sensitivity) > 0 else 0.0
    return sensitivity, precision, f1


def evaluate(slt_variants, truth_snv, truth_indel, regions_by_chrom):
    """Evaluate SLT classifications against truth sets.

    For each cumulative tier threshold (A, A+B, A+B+C, A+B+C+D):
      - TP: variant in SLT set AND in truth set AND in evaluation regions
      - FP: variant in SLT set AND NOT in truth set AND in evaluation regions
      - FN: variant in truth set AND NOT in SLT set AND in evaluation regions
    """
    truth_all = truth_snv | truth_indel
    tier_order = ["SLT-A", "SLT-B", "SLT-C", "SLT-D"]

    # Filter SLT variants to evaluation regions
    slt_in_regions = [v for v in slt_variants
                      if variant_in_regions(v["chrom"], v["pos"], regions_by_chrom)]

    # Filter truth to evaluation regions (should already be, but verify)
    truth_in_regions = {t for t in truth_all
                        if variant_in_regions(t[0], t[1], regions_by_chrom)}

    results = []

    # Cumulative evaluation: include tiers up to threshold
    included_tiers = set()
    for tier in tier_order:
        included_tiers.add(tier)

        # Variants called somatic at this threshold
        called = {v["key"] for v in slt_in_regions if v["tier"] in included_tiers}

        tp = len(called & truth_in_regions)
        fp = len(called - truth_in_regions)
        fn = len(truth_in_regions - called)

        sens, prec, f1 = compute_metrics(tp, fp, fn)

        results.append({
            "threshold": f"≥{tier}",
            "tiers_included": "+".join(sorted(included_tiers)),
            "n_called": len(called),
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "sensitivity": sens,
            "precision": prec,
            "f1": f1,
        })

    # Also evaluate per-tier (exclusive)
    for tier in tier_order:
        tier_variants = {v["key"] for v in slt_in_regions if v["tier"] == tier}
        tp = len(tier_variants & truth_in_regions)
        fp = len(tier_variants - truth_in_regions)

        results.append({
            "threshold": f"{tier}_only",
            "tiers_included": tier,
            "n_called": len(tier_variants),
            "tp": tp,
            "fp": fp,
            "fn": "N/A",
            "sensitivity": "N/A",
            "precision": tp / (tp + fp) if (tp + fp) > 0 else 0.0,
            "f1": "N/A",
        })

    # SNV-only evaluation
    slt_snv = [v for v in slt_in_regions
               if len(v["ref"]) == 1 and len(v["alt"]) == 1]
    truth_snv_in_regions = {t for t in truth_snv
                            if variant_in_regions(t[0], t[1], regions_by_chrom)}

    included_tiers = set()
    for tier in tier_order:
        included_tiers.add(tier)
        called = {v["key"] for v in slt_snv if v["tier"] in included_tiers}
        tp = len(called & truth_snv_in_regions)
        fp = len(called - truth_snv_in_regions)
        fn = len(truth_snv_in_regions - called)
        sens, prec, f1 = compute_metrics(tp, fp, fn)

        results.append({
            "threshold": f"SNV_≥{tier}",
            "tiers_included": "+".join(sorted(included_tiers)),
            "n_called": len(called),
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "sensitivity": sens,
            "precision": prec,
            "f1": f1,
        })

    return results, len(slt_in_regions), len(truth_in_regions), len(truth_snv_in_regions)


def write_per_variant(slt_variants, truth_all, regions_by_chrom, output_path):
    """Write per-variant evaluation table."""
    fields = ["chrom", "pos", "ref", "alt", "tier", "posterior", "gene",
              "consequence", "filter", "in_truth", "in_eval_region"]

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()

        for v in slt_variants:
            in_region = variant_in_regions(v["chrom"], v["pos"], regions_by_chrom)
            in_truth = v["key"] in truth_all

            writer.writerow({
                "chrom": v["chrom"],
                "pos": v["pos"],
                "ref": v["ref"],
                "alt": v["alt"],
                "tier": v["tier"],
                "posterior": v["posterior"],
                "gene": v["gene"],
                "consequence": v["consequence"],
                "filter": v["filter"],
                "in_truth": in_truth,
                "in_eval_region": in_region,
            })


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate SLT classifications against SEQC2 truth set"
    )
    parser.add_argument("--slt-input", required=True, help="SLT classified TSV")
    parser.add_argument("--truth-snv", required=True, help="Truth SNV VCF (gzipped)")
    parser.add_argument("--truth-indel", required=True, help="Truth indel VCF (gzipped)")
    parser.add_argument("--eval-bed", required=True, help="Evaluation regions BED")
    parser.add_argument("--output", required=True, help="Per-variant evaluation TSV")
    parser.add_argument("--summary", required=True, help="Summary metrics text file")
    args = parser.parse_args()

    print("Loading truth sets...")
    truth_snv = parse_vcf_truth(args.truth_snv)
    truth_indel = parse_vcf_truth(args.truth_indel)
    print(f"  Truth SNVs: {len(truth_snv)}")
    print(f"  Truth indels: {len(truth_indel)}")

    print("Loading evaluation regions...")
    regions = parse_bed_regions(args.eval_bed)
    regions_by_chrom = build_region_index(regions)
    print(f"  Evaluation regions: {len(regions)}")

    print("Loading SLT classifications...")
    slt_variants = parse_slt_variants(args.slt_input)
    print(f"  SLT variants: {len(slt_variants)}")

    # Tier distribution
    tier_counts = defaultdict(int)
    for v in slt_variants:
        tier_counts[v["tier"]] += 1
    for tier in ["SLT-A", "SLT-B", "SLT-C", "SLT-D"]:
        print(f"    {tier}: {tier_counts.get(tier, 0)}")

    print("Evaluating...")
    truth_all = truth_snv | truth_indel
    results, n_slt_eval, n_truth_eval, n_truth_snv_eval = evaluate(
        slt_variants, truth_snv, truth_indel, regions_by_chrom
    )

    # Write per-variant table
    print("Writing per-variant evaluation...")
    write_per_variant(slt_variants, truth_all, regions_by_chrom, args.output)

    # Write summary
    with open(args.summary, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("SEQC2 HCC1395 Benchmark — SLT Evaluation Summary\n")
        f.write("=" * 70 + "\n\n")

        f.write(f"Truth SNVs (total): {len(truth_snv)}\n")
        f.write(f"Truth indels (total): {len(truth_indel)}\n")
        f.write(f"Truth in eval regions: {n_truth_eval}\n")
        f.write(f"Truth SNVs in eval regions: {n_truth_snv_eval}\n")
        f.write(f"SLT variants in eval regions: {n_slt_eval}\n\n")

        f.write("-" * 70 + "\n")
        f.write(f"{'Threshold':<20} {'Called':>8} {'TP':>8} {'FP':>8} {'FN':>8} "
                f"{'Sens':>8} {'Prec':>8} {'F1':>8}\n")
        f.write("-" * 70 + "\n")

        for r in results:
            sens = f"{r['sensitivity']:.4f}" if isinstance(r['sensitivity'], float) else r['sensitivity']
            prec = f"{r['precision']:.4f}" if isinstance(r['precision'], float) else r['precision']
            f1 = f"{r['f1']:.4f}" if isinstance(r['f1'], float) else r['f1']
            f.write(f"{r['threshold']:<20} {r['n_called']:>8} {r['tp']:>8} "
                    f"{str(r['fp']):>8} {str(r['fn']):>8} "
                    f"{sens:>8} {prec:>8} {f1:>8}\n")

        f.write("-" * 70 + "\n")

    # Print summary to stdout too
    with open(args.summary) as f:
        print(f.read())

    print(f"\nPer-variant output: {args.output}")
    print(f"Summary: {args.summary}")


if __name__ == "__main__":
    main()
