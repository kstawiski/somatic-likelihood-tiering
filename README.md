# Somatic Likelihood Tiering (SLT)

An interpretable, open-source post-calling triage framework for somatic variant classification in tumor-only whole-exome sequencing (WES).

## Overview

SLT classifies variants from tumor-only WES into four confidence tiers (SLT-A through SLT-D) using:

- **Four complementary evidence layers**: population allele frequency (POPAF), gnomAD annotation, germline quality (GERMQ), and COSMIC recurrence
- **PureCN posterior somatic probabilities**: Bayesian posterior from copy number-aware classification
- **Integrated CHIP detection**: Clonal hematopoiesis of indeterminate potential flagging with curated gene lists

All thresholds are deterministic, convention-grounded, frozen, and fully auditable.

## Key Features

- **Zero dependencies**: Core classifier uses only the Python standard library (Python 3.7+)
- **Interpretable**: Every tier assignment is traceable to specific evidence layers
- **Graduated confidence**: Four tiers for different downstream applications (clinical reporting, discovery, screening)
- **CHIP-aware**: Integrated clonal hematopoiesis detection prevents false somatic calls from blood-derived variants
- **Caller-agnostic**: Works with Mutect2, VarDict, or any TSV/MAF-producing pipeline
- **Annotation-only mode**: Operates with public databases only (gnomAD + COSMIC) when PureCN or BAM-level annotations are unavailable

## Installation

### Standalone (no installation needed)

```bash
# Download and run directly - no pip install required
python3 slt_classify.py --input variants.tsv --output classified.tsv
```

### pip install

```bash
pip install .
# Then use the CLI:
slt-classify --input variants.tsv --output classified.tsv
```

### Development

```bash
pip install -e ".[test]"
pytest tests/ -v
```

## Quick Start

### From a MAF file (most common)

```bash
# 1. Convert MAF to SLT input format
python maf_to_slt_input.py --input somatic.maf --output slt_input.tsv

# 2. Run SLT (annotation-only mode if no PureCN)
python slt_classify.py --input slt_input.tsv --output classified.tsv --annotation-only

# Or with PureCN posteriors for full mode:
python maf_to_slt_input.py --input somatic.maf --purecn purecn_calls.csv --output slt_input.tsv
python slt_classify.py --input slt_input.tsv --output classified.tsv
```

### From a VCF file

```bash
# 1. Convert annotated VCF to SLT input format
python vcf_to_slt_input.py --input annotated.vcf.gz --output slt_input.tsv

# 2. Run SLT
python slt_classify.py --input slt_input.tsv --output classified.tsv
```

### Direct usage (pre-formatted TSV)

```bash
# Full mode (Mutect2 + PureCN)
python slt_classify.py --input annotated_variants.tsv --output slt_output.tsv

# Annotation-only mode (MAF-level, no PureCN)
python slt_classify.py --input variants.maf --output slt_output.tsv --annotation-only
```

### Python API

```python
from slt_classify import classify_variant

variant = {
    "POSTERIOR.SOMATIC": "0.95",
    "POPAF": "40",
    "gnomAD_AF": "0.00001",
    "GERMQ": "50",
    "COSMIC_CONFIRMED_SOMATIC": "20",
    "Hugo_Symbol": "FGFR3",
    "t_alt_freq": "0.30",
    "HGVSp_Short": "S249C",
    "Variant_Classification": "Missense_Mutation",
}

result = classify_variant(variant)
print(result["slt_tier"])           # SLT-A
print(result["slt_evidence_level"]) # high
print(result["slt_chip_status"])    # no_chip
```

## Preparing Input from Common Formats

SLT accepts a simple tab-separated input, but most variant callers produce VCF or MAF files. Two helper scripts are provided to convert these formats.

### From MAF files (`maf_to_slt_input.py`)

Converts standard MAF (Mutation Annotation Format) files to SLT input. Handles column naming conventions from Funcotator, Oncotator, maf2maf, and cBioPortal.

```bash
# Basic conversion (annotation-only mode — no PureCN)
python maf_to_slt_input.py --input somatic.maf --output slt_input.tsv

# With PureCN posteriors for full-mode SLT
python maf_to_slt_input.py --input somatic.maf --purecn purecn_calls.csv --output slt_input.tsv
```

**Features:**
- Auto-detects column names (Hugo_Symbol vs GENE vs gene, etc.)
- Computes VAF from t_ref_count/t_alt_count if no direct VAF column
- Merges PureCN posteriors, POPAF, and GERMQ by genomic coordinates
- Extracts COSMIC counts from existing annotation fields

### From VCF files (`vcf_to_slt_input.py`)

Converts annotated VCF files (Funcotator or VEP) to SLT input. The VCF must be annotated before conversion — see [Reference Files](#reference-files-for-annotation) below.

```bash
# From Funcotator-annotated VCF
python vcf_to_slt_input.py --input funcotator_annotated.vcf --output slt_input.tsv

# From VEP-annotated VCF
python vcf_to_slt_input.py --input vep_annotated.vcf.gz --output slt_input.tsv

# With PureCN posteriors
python vcf_to_slt_input.py --input annotated.vcf --purecn purecn_calls.csv --output slt_input.tsv

# Specify sample (for multi-sample VCFs)
python vcf_to_slt_input.py --input annotated.vcf --sample TUMOR --output slt_input.tsv
```

**Features:**
- Parses Funcotator FUNCOTATION and VEP CSQ INFO fields
- Extracts Mutect2 POPAF and GERMQ from INFO
- Computes VAF from FORMAT/AD fields
- Handles gzipped VCFs (.vcf.gz)

## Reference Files for Annotation

SLT itself has no reference file dependencies — it operates on pre-annotated TSV input. However, the upstream annotation pipeline requires several reference databases. Below is the recommended annotation workflow.

### Required for Full Mode

| Resource | Version | Purpose | Download |
|----------|---------|---------|----------|
| **Reference genome** | GRCh38 / hg38 | Alignment, variant calling | [GATK bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811) |
| **Funcotator data sources** | v1.8 (hg38) | Gene annotation, gnomAD AF, COSMIC | `gatk FuncotatorDataSourceDownloader --germline --validate-integrity --extract-after-download` |
| **gnomAD** | v4.1 exomes | Population allele frequencies (Layer 2) | [gnomAD downloads](https://gnomad.broadinstitute.org/downloads) |
| **COSMIC** | v103+ | Somatic recurrence (Layer 4) | [COSMIC](https://cancer.sanger.ac.uk/cosmic/download) (registration required) |
| **PureCN panel of normals** | Project-specific | Posterior somatic probability | Built from matched normals or unmatched panel |

### Deployment Guidelines

| Setting | Recommended Mode | Maximum Tier | Notes |
|---------|-----------------|--------------|-------|
| Process-matched NormalDB available | Full SLT | SLT-A | Best performance; requires BAM files + matched normals |
| No NormalDB / no BAM access | Annotation-only | SLT-C | Uses only gnomAD + COSMIC; `--annotation-only` flag |
| Matched normal tissue available | Paired T/N analysis | N/A | Preferred over any tumor-only approach |

### Recommended Annotation Pipeline

```bash
# 1. Variant calling (Mutect2 tumor-only)
gatk Mutect2 \
    -R reference.fa \
    -I tumor.bam \
    --germline-resource af-only-gnomad.hg38.vcf.gz \
    --panel-of-normals pon.vcf.gz \
    -O raw.vcf.gz

# 2. Filter variants
gatk FilterMutectCalls \
    -R reference.fa \
    -V raw.vcf.gz \
    --contamination-table contamination.table \
    -O filtered.vcf.gz

# 3. Annotate with Funcotator (provides gene, consequence, gnomAD, COSMIC)
gatk Funcotator \
    -R reference.fa \
    -V filtered.vcf.gz \
    --ref-version hg38 \
    --data-sources-path funcotator_dataSources.v1.8.hg38.20230908s \
    -O annotated.vcf \
    --output-file-format VCF

# 4. Run PureCN (provides posterior somatic probability)
# See PureCN documentation: https://bioconductor.org/packages/PureCN/
Rscript PureCN.R \
    --sampleid SAMPLE \
    --tumor filtered.vcf.gz \
    --normaldb normalDB.rds \
    --intervals targets_intervals.txt \
    --genome hg38

# 5. Convert to SLT input and classify
python vcf_to_slt_input.py \
    --input annotated.vcf \
    --purecn PureCN_output/SAMPLE.csv \
    --output slt_input.tsv

python slt_classify.py --input slt_input.tsv --output slt_classified.tsv
```

### Annotation-Only Mode (Minimal Requirements)

For MAF-level reanalysis without BAM files, only gnomAD and COSMIC annotations are needed:

| Resource | Purpose | SLT Layer |
|----------|---------|-----------|
| **gnomAD AF** | Population frequency filtering | Layer 2 |
| **COSMIC** | Somatic recurrence database | Layer 4 |

Layers 1 (POPAF) and 3 (GERMQ) require Mutect2 BAM-level output and are disabled in annotation-only mode. PureCN posterior is unavailable. Maximum achievable tier is SLT-C.

### Alternative Annotation with VEP

```bash
# Ensembl VEP (alternative to Funcotator)
vep --input_file filtered.vcf.gz \
    --output_file annotated.vcf \
    --format vcf --vcf \
    --assembly GRCh38 \
    --cache --dir_cache /path/to/vep_cache \
    --fasta reference.fa \
    --everything --pick \
    --plugin COSMIC,/path/to/CosmicCodingMuts.vcf.gz \
    --custom gnomAD_exomes.vcf.gz,gnomADe,vcf,exact,0,AF
```

## Input Format

Tab-separated file with the following columns (flexible naming):

| Field | Alternative Names | Description |
|-------|-------------------|-------------|
| `POSTERIOR.SOMATIC` | `POSTERIOR_SOMATIC` | PureCN posterior somatic probability |
| `POPAF` | | Mutect2 negative log10 population allele frequency |
| `gnomAD_AF` | `gnomAD_exome_AF` | gnomAD exome allele frequency |
| `GERMQ` | | Mutect2 germline quality score |
| `COSMIC_CONFIRMED_SOMATIC` | | COSMIC confirmed somatic count |
| `COSMIC_HOTSPOT` | | COSMIC hotspot count |
| `COSMIC_SAMPLE` | `COSMIC_TOTAL_OCC` | COSMIC sample count |
| `CGC_GENE` | | Cancer Gene Census membership (TRUE/FALSE) |
| `Hugo_Symbol` | `GENE`, `gene` | Gene symbol |
| `t_alt_freq` | `AF`, `VAF` | Variant allele frequency |
| `HGVSp_Short` | `Protein_Change`, `AAChange` | Protein change |
| `Variant_Classification` | `Consequence` | Variant consequence |

Missing fields are handled gracefully (treated as absent evidence).

## Output

The classifier appends 10 columns to each input row:

| Column | Type | Description |
|--------|------|-------------|
| `slt_layer1_popaf` | bool | POPAF >= 5.0 |
| `slt_layer2_gnomad` | bool | gnomAD AF < 0.001 or absent |
| `slt_layer3_germq` | bool | GERMQ >= 30 |
| `slt_layer4_cosmic` | bool | COSMIC composite criterion |
| `slt_n_somatic_layers` | int | Number of passing layers (0-4) |
| `slt_evidence_level` | str | high / medium / low |
| `slt_chip_status` | str | chip_likely / chip_possible / no_chip |
| `slt_posterior_used` | float | PureCN posterior (or "NA") |
| `slt_tier` | str | SLT-A / SLT-B / SLT-C / SLT-D |
| `slt_mode` | str | full / annotation_only |

## Tier Definitions

### Classification Cascade

Conditions are evaluated in order; the first match determines the tier:

| Tier | Meaning | Conditions |
|------|---------|------------|
| **SLT-A** | High confidence somatic | Posterior >= 0.8 AND evidence in {high, medium} AND not chip_likely |
| **SLT-B** | Likely somatic | Posterior >= 0.5, OR (high evidence AND >= 3 layers); chip_likely blocked |
| **SLT-C** | Possible somatic | Posterior >= 0.2, OR evidence in {high, medium} |
| **SLT-D** | Unlikely somatic | All remaining variants |

### Clinical Workflow

| Tier | Recommended Action |
|------|-------------------|
| **SLT-A** | Report as somatic; include in TMB; prioritize for targeted therapy matching |
| **SLT-B** | Report with annotation; manual review recommended; consider orthogonal confirmation for treatment-critical variants |
| **SLT-C** | Include in discovery analyses; flag for extended review if clinically relevant gene |
| **SLT-D** | Deprioritize; do not include in TMB |
| **CHIP-likely** | Separate reporting; flag for hematology follow-up if clinically indicated |

### Evidence Levels

| Level | Criteria |
|-------|----------|
| **High** | >= 3 somatic-supporting layers AND COSMIC layer passes |
| **Medium** | >= 2 somatic-supporting layers |
| **Low** | < 2 somatic-supporting layers |

### CHIP Classification

Variants in CHIP-associated genes are evaluated for clonal hematopoiesis:

- **Tier 1 genes** (13 canonical CHIP drivers): DNMT3A, TET2, ASXL1, PPM1D, JAK2, SF3B1, SRSF2, U2AF1, IDH1, IDH2, ZBTB33, GNB1, CBL
- **Tier 2 genes** (28 extended): TP53, KRAS, NRAS, FLT3, KIT, NPM1, and others

`chip_likely` status blocks both SLT-A and SLT-B assignment, preventing CHIP variants from reaching the most actionable tiers. `chip_likely` variants are downgraded to SLT-C (maximum) or SLT-D. Tier 2 genes (including TP53, KRAS, NRAS) receive `chip_possible` annotation but are NOT blocked from SLT-A/B.

### Annotation-Only Mode

When PureCN is unavailable (e.g., MAF-level reanalysis without BAM files):

- Layers 1 (POPAF) and 3 (GERMQ) are disabled
- PureCN posterior is set to null
- Maximum achievable tier is SLT-C
- SLT-A and SLT-B are unreachable

Use the `--annotation-only` flag (or `--degraded` for backward compatibility).

## Thresholds

All thresholds are convention-grounded, frozen, and deterministic:

| Parameter | Value | Source |
|-----------|-------|--------|
| Posterior SLT-A gate | >= 0.8 | PureCN recommended (Riester et al., 2016) |
| Posterior SLT-B gate | >= 0.5 | PureCN recommended |
| Posterior SLT-C gate | >= 0.2 | PureCN recommended |
| POPAF threshold | >= 5.0 | Mutect2 standard (AF <= 1e-5) |
| gnomAD AF threshold | < 0.001 | ACMG/AMP BA1/BS1 aligned |
| GERMQ threshold | >= 30 | Mutect2 standard (<=0.1% germline probability) |
| COSMIC confirmed min | >= 5 | Conservative recurrence threshold |
| COSMIC hotspot min | >= 10 | Conservative recurrence threshold |
| CGC + COSMIC sample min | >= 2 | Conservative recurrence threshold |

## Benchmark Performance

### SEQC2 HCC1395 (primary benchmark)

Validated on the SEQC2 HCC1395 breast cancer truth set (455 variants in evaluation regions), processed in tumor-only mode:

| Threshold | Called | TP | FP | Sensitivity | PPV | F1 | NNR |
|-----------|-------|----|----|-------------|-----|-----|-----|
| >= SLT-A | 105 | 82 | 23 | 18.0% | 78.1% | 0.293 | 1.28 |
| >= SLT-B | 160 | 102 | 58 | 22.4% | 63.7% | 0.332 | 1.57 |
| >= SLT-C | 2,423 | 422 | 2,001 | 92.7% | 17.4% | 0.293 | 5.74 |

PureCN posterior AUROC: 0.775 (95% CI: 0.732-0.817). 100% CGC driver retention at SLT-C. Evidence layers rescued 77.3% of true positives lacking PureCN support.

### External validation (BostonGene cell lines)

| Sample | Tumor type | Truth (n) | SLT >=C Sensitivity | SLT >=C PPV | SLT-attributable FN |
|--------|-----------|-----------|---------------------|-------------|---------------------|
| COLO829 | Melanoma | 357 | 87.1% | 5.1% | 2 (0.6%) |
| NCI-H1770 | NSCLC | 1,042 | 81.9% | 7.9% | 23 (2.2%) |

### Clinical validation (HdM-BLCA-1 bladder cancer)

Validated on 22 FFPE metastatic urothelial carcinoma patients (Boll et al. 2023, *Sci Rep*; EGA: EGAS00001007086). Reference standard: matched tumor-normal Mutect2 PASS variants (TLOD >= 20, VAF >= 5%).

| Threshold | Sensitivity [95% CI] | PPV | NNR |
|-----------|---------------------|-----|-----|
| >= SLT-A | 0.09% [0.01-0.24%] | 0.07% | 1,411 |
| >= SLT-B | 0.50% [0.34-0.72%] | 0.30% | 336 |
| >= SLT-C | 60.0% [51.7-69.6%] | 0.69% | 145 |
| All tiers | 94.4% [88.7-98.6%] | 0.63% | 158 |

All-tier sensitivity (94.4%) matches the SEQC2 cell-line benchmark (94.5%), confirming that the evidence-layer architecture generalizes to clinical FFPE samples. SLT-A/B precision tiers showed limited PureCN posterior coverage on FFPE (<5% of variants), establishing a deployment boundary for FFPE workflows. Per-patient SLT-C sensitivity: median 69.7% (IQR 58.7-74.6%).

## Tests

```bash
# Install test dependencies
pip install pytest

# Run tests
pytest tests/ -v
```

The test suite covers all evidence layers, CHIP classification, tier assignment cascade, annotation-only mode, boundary conditions, and integration tests.

## Repository Structure

```
somatic-likelihood-tiering/
├── slt_classify.py          # Core classifier (standalone, no dependencies)
├── maf_to_slt_input.py      # MAF → SLT input converter
├── vcf_to_slt_input.py      # VCF → SLT input converter (Funcotator/VEP)
├── tests/
│   └── test_slt_classify.py # Unit test suite (39+ tests)
├── examples/
│   ├── example_input.tsv    # Example input
│   └── example_output.tsv   # Example output
├── benchmark/               # Benchmark evaluation scripts
│   └── README.md
├── pyproject.toml           # Python packaging
├── LICENSE                  # MIT License
└── README.md                # This file
```

## Citation

If you use SLT in your research, please cite:

> Stawiski K, Kamran SC, De Carvalho FLF, Mouw KW. Somatic Likelihood Tiering (SLT): an interpretable post-calling triage framework for tumor-only whole-exome sequencing. *BMC Bioinformatics* (under review), 2026.

## License

MIT License. See [LICENSE](LICENSE) for details.
