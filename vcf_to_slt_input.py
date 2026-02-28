#!/usr/bin/env python3
"""Convert an annotated VCF to SLT input format.

Reads a VCF file annotated with Funcotator or VEP and produces a
tab-separated file compatible with slt_classify.py.

Supports:
- Funcotator-annotated VCFs (FUNCOTATION INFO field)
- VEP-annotated VCFs (CSQ INFO field)
- Mutect2 VCFs (POPAF, GERMQ from FilterMutectCalls)
- PureCN-annotated VCFs (optional POSTERIOR.SOMATIC)

No external dependencies required (Python 3.7+ standard library only).

Usage:
    python vcf_to_slt_input.py --input annotated.vcf --output slt_input.tsv
    python vcf_to_slt_input.py --input annotated.vcf.gz --output slt_input.tsv --purecn purecn_calls.csv

Reference files for annotation (before running this script):
    # Funcotator annotation (recommended)
    gatk Funcotator \\
        -V mutect2_filtered.vcf.gz \\
        -R reference.fa \\
        --ref-version hg38 \\
        --data-sources-path funcotator_dataSources.v1.8.hg38.20230908s \\
        -O annotated.vcf \\
        --output-file-format VCF

    # VEP annotation (alternative)
    vep --input_file mutect2_filtered.vcf.gz \\
        --output_file annotated.vcf \\
        --format vcf --vcf \\
        --assembly GRCh38 \\
        --cache --dir_cache /path/to/vep_cache \\
        --fasta reference.fa \\
        --everything --pick
"""

import argparse
import csv
import gzip
import re
import sys
from typing import Dict, List, Optional, Tuple


SLT_COLUMNS = [
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Variant_Classification",
    "HGVSp_Short",
    "t_alt_freq",
    "POPAF",
    "GERMQ",
    "gnomAD_AF",
    "COSMIC_CONFIRMED_SOMATIC",
    "COSMIC_HOTSPOT",
    "COSMIC_SAMPLE",
    "CGC_GENE",
    "POSTERIOR.SOMATIC",
    "Tumor_Sample_Barcode",
]


def open_vcf(path: str):
    """Open a VCF file (handles .gz)."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_info(info_str: str) -> dict:
    """Parse VCF INFO field into dict."""
    info = {}
    for item in info_str.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            info[key] = val
        else:
            info[item] = True
    return info


def parse_funcotation(funcotation_str: str, funcotation_keys: List[str]) -> dict:
    """Parse Funcotator FUNCOTATION field."""
    # FUNCOTATION=[val1|val2|val3|...]
    clean = funcotation_str.strip("[]")
    # Multiple transcripts separated by ],[ — take first
    if "],[" in clean:
        clean = clean.split("],[")[0]
    values = clean.split("|")
    result = {}
    for i, key in enumerate(funcotation_keys):
        if i < len(values):
            val = values[i]
            result[key] = val if val not in ("", "__UNKNOWN__") else ""
    return result


def parse_csq(csq_str: str, csq_keys: List[str]) -> dict:
    """Parse VEP CSQ field (take first/picked transcript)."""
    transcripts = csq_str.split(",")
    # Prefer PICK=1 transcript if available
    best = transcripts[0]
    for t in transcripts:
        vals = t.split("|")
        # Check for PICK flag
        if "PICK" in csq_keys:
            pick_idx = csq_keys.index("PICK")
            if pick_idx < len(vals) and vals[pick_idx] == "1":
                best = t
                break

    values = best.split("|")
    result = {}
    for i, key in enumerate(csq_keys):
        if i < len(values):
            result[key] = values[i]
    return result


def extract_format_field(format_str: str, sample_str: str, field: str) -> Optional[str]:
    """Extract a field from VCF FORMAT/SAMPLE columns."""
    keys = format_str.split(":")
    vals = sample_str.split(":")
    if field in keys:
        idx = keys.index(field)
        if idx < len(vals):
            val = vals[idx]
            return val if val != "." else None
    return None


def compute_vaf_from_ad(format_str: str, sample_str: str) -> Optional[str]:
    """Compute VAF from AD (allelic depth) field."""
    ad = extract_format_field(format_str, sample_str, "AD")
    if ad:
        counts = ad.split(",")
        if len(counts) >= 2:
            try:
                ref = int(counts[0])
                alt = int(counts[1])
                total = ref + alt
                if total > 0:
                    return f"{alt / total:.6f}"
            except ValueError:
                pass
    # Try AF field
    af = extract_format_field(format_str, sample_str, "AF")
    if af:
        return af
    return None


def load_purecn(purecn_path: str) -> Dict[str, dict]:
    """Load PureCN output and index by chr:pos:ref:alt."""
    purecn_data = {}
    with open(purecn_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            chrom = row.get("chr", "")
            pos = row.get("start", "")
            ref = row.get("ref", "")
            alt = row.get("alt", "")
            key = f"{chrom}:{pos}:{ref}:{alt}"
            purecn_data[key] = {
                "POSTERIOR.SOMATIC": row.get("POSTERIOR.SOMATIC", ""),
                "POPAF": row.get("POPAF", ""),
                "GERMQ": row.get("GERMQ", ""),
            }
    return purecn_data


def vep_consequence_to_maf(consequence: str) -> str:
    """Map VEP SO terms to MAF Variant_Classification (best effort)."""
    # Only map the primary term (first in & or ,-separated list)
    mapping = {
        "missense_variant": "Missense_Mutation",
        "nonsense": "Nonsense_Mutation",
        "stop_gained": "Nonsense_Mutation",
        "frameshift_variant": "Frame_Shift_Del",
        "splice_donor_variant": "Splice_Site",
        "splice_acceptor_variant": "Splice_Site",
        "splice_region_variant": "Splice_Region",
        "start_lost": "Translation_Start_Site",
        "stop_lost": "Nonstop_Mutation",
        "inframe_insertion": "In_Frame_Ins",
        "inframe_deletion": "In_Frame_Del",
        "synonymous_variant": "Silent",
        "5_prime_UTR_variant": "5'UTR",
        "3_prime_UTR_variant": "3'UTR",
        "intron_variant": "Intron",
        "upstream_gene_variant": "5'Flank",
        "downstream_gene_variant": "3'Flank",
        "intergenic_variant": "IGR",
    }
    return mapping.get(consequence, consequence)


def main():
    parser = argparse.ArgumentParser(
        description="Convert annotated VCF to SLT input format",
        epilog="Supports Funcotator and VEP annotations",
    )
    parser.add_argument("--input", "-i", required=True, help="Input VCF (or .vcf.gz)")
    parser.add_argument("--output", "-o", required=True, help="Output SLT-ready TSV")
    parser.add_argument("--purecn", help="Optional PureCN output CSV to merge posteriors")
    parser.add_argument("--sample", help="Sample name (default: first sample in VCF)")
    args = parser.parse_args()

    purecn_data = {}
    if args.purecn:
        print(f"Loading PureCN data from {args.purecn}...", file=sys.stderr)
        purecn_data = load_purecn(args.purecn)
        print(f"  {len(purecn_data)} variants loaded", file=sys.stderr)

    funcotation_keys = []
    csq_keys = []
    has_funcotator = False
    has_vep = False
    sample_name = args.sample or ""
    sample_idx = 0

    with open_vcf(args.input) as fin, \
         open(args.output, "w", newline="") as fout:

        writer = csv.DictWriter(fout, fieldnames=SLT_COLUMNS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()

        n = 0
        n_purecn = 0

        for line in fin:
            line = line.strip()

            # Parse header
            if line.startswith("##"):
                # Check for Funcotator keys
                m = re.match(r'##INFO=<ID=FUNCOTATION,.*Description="Funcotation fields are: (.+)">', line)
                if m:
                    funcotation_keys = m.group(1).split("|")
                    has_funcotator = True
                    print(f"  Funcotator annotation detected ({len(funcotation_keys)} fields)", file=sys.stderr)

                # Check for VEP CSQ keys
                m = re.match(r'##INFO=<ID=CSQ,.*Format: (.+)">', line)
                if m:
                    csq_keys = m.group(1).split("|")
                    has_vep = True
                    print(f"  VEP annotation detected ({len(csq_keys)} fields)", file=sys.stderr)
                continue

            if line.startswith("#CHROM"):
                cols = line.split("\t")
                if len(cols) > 9:
                    if args.sample:
                        if args.sample in cols[9:]:
                            sample_idx = cols.index(args.sample) - 9
                            sample_name = args.sample
                        else:
                            print(f"WARNING: Sample '{args.sample}' not found in VCF; using first sample",
                                  file=sys.stderr)
                            sample_name = cols[9]
                    else:
                        sample_name = cols[9]
                    print(f"  Sample: {sample_name}", file=sys.stderr)
                continue

            # Data line
            fields = line.split("\t")
            if len(fields) < 8:
                continue

            chrom, pos, _, ref, alt, _, filt, info_str = fields[:8]
            format_str = fields[8] if len(fields) > 8 else ""
            sample_str = fields[9 + sample_idx] if len(fields) > 9 + sample_idx else ""

            # Handle multi-allelic: take first ALT
            if "," in alt:
                alt = alt.split(",")[0]

            info = parse_info(info_str)

            out = {col: "" for col in SLT_COLUMNS}
            out["Chromosome"] = chrom.replace("chr", "")
            out["Start_Position"] = pos
            out["Reference_Allele"] = ref
            out["Tumor_Seq_Allele2"] = alt
            out["Tumor_Sample_Barcode"] = sample_name

            # Extract Mutect2 INFO fields
            if "POPAF" in info:
                out["POPAF"] = str(info["POPAF"])
            if "GERMQ" in info:
                out["GERMQ"] = str(info["GERMQ"])

            # Extract annotation
            if has_funcotator and "FUNCOTATION" in info:
                func = parse_funcotation(str(info["FUNCOTATION"]), funcotation_keys)
                out["Hugo_Symbol"] = func.get("Gencode_34_hugoSymbol", "") or func.get("Gencode_19_hugoSymbol", "")
                out["Variant_Classification"] = func.get("Gencode_34_variantClassification", "") or \
                                                 func.get("Gencode_19_variantClassification", "")
                out["HGVSp_Short"] = func.get("Gencode_34_proteinChange", "") or \
                                      func.get("Gencode_19_proteinChange", "")
                # gnomAD from Funcotator
                gnomad_val = func.get("gnomAD_exome_AF", "") or func.get("gnomAD_AF", "")
                if gnomad_val:
                    out["gnomAD_AF"] = gnomad_val
                # COSMIC
                cosmic_val = func.get("COSMIC_overlapping_mutations", "")
                if cosmic_val:
                    # Count occurrences
                    out["COSMIC_SAMPLE"] = str(cosmic_val.count("|") + 1) if cosmic_val else ""

            elif has_vep and "CSQ" in info:
                csq = parse_csq(str(info["CSQ"]), csq_keys)
                out["Hugo_Symbol"] = csq.get("SYMBOL", "")
                consequence = csq.get("Consequence", "")
                out["Variant_Classification"] = consequence
                hgvsp = csq.get("HGVSp", "")
                if hgvsp and ":" in hgvsp:
                    hgvsp = hgvsp.split(":")[-1]
                out["HGVSp_Short"] = hgvsp
                gnomad_val = csq.get("gnomADe_AF", "") or csq.get("gnomAD_AF", "")
                if gnomad_val:
                    out["gnomAD_AF"] = gnomad_val
                # COSMIC from VEP
                existing = csq.get("Existing_variation", "")
                if existing and "COSV" in existing:
                    cosm_ids = [v for v in existing.split("&") if v.startswith("COSV") or v.startswith("COSM")]
                    if cosm_ids:
                        out["COSMIC_SAMPLE"] = str(len(cosm_ids))

            # Compute VAF from FORMAT/SAMPLE
            if format_str and sample_str:
                vaf = compute_vaf_from_ad(format_str, sample_str)
                if vaf:
                    out["t_alt_freq"] = vaf

            # Merge PureCN data
            if purecn_data:
                key = f"{chrom}:{pos}:{ref}:{alt}"
                # Also try with chr prefix stripped/added
                key_alt = f"chr{chrom}:{pos}:{ref}:{alt}" if not chrom.startswith("chr") else \
                          f"{chrom.replace('chr','')}:{pos}:{ref}:{alt}"
                pc = purecn_data.get(key) or purecn_data.get(key_alt)
                if pc:
                    for field in ["POSTERIOR.SOMATIC", "POPAF", "GERMQ"]:
                        if pc.get(field) and not out.get(field):
                            out[field] = pc[field]
                    n_purecn += 1

            writer.writerow(out)
            n += 1

    print(f"\nConverted {n} variants to SLT input format", file=sys.stderr)
    if purecn_data:
        print(f"  PureCN matched: {n_purecn}/{n}", file=sys.stderr)
    print(f"  Output: {args.output}", file=sys.stderr)

    mode_hint = ""
    if not purecn_data:
        if has_funcotator:
            print("\n  Tip: Funcotator VCF detected. For full SLT mode, also provide --purecn.",
                  file=sys.stderr)
            mode_hint = " --degraded"
        elif has_vep:
            print("\n  Tip: VEP VCF detected. For full SLT mode, also provide --purecn.",
                  file=sys.stderr)
            mode_hint = " --degraded"

    print(f"\nNext step: python slt_classify.py --input {args.output} --output classified.tsv{mode_hint}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
