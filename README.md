# Somatic Likelihood Tiering (SLT)

An interpretable, multi-evidence framework for somatic variant classification in tumor-only whole-exome sequencing (WES).

## Overview

SLT classifies variants from tumor-only WES into four confidence tiers (SLT-A through SLT-D) using:

- **Four complementary evidence layers**: population allele frequency (POPAF), gnomAD annotation, germline quality (GERMQ), and COSMIC recurrence
- **PureCN posterior somatic probabilities**: Bayesian posterior from copy number-aware classification
- **Integrated CHIP detection**: Clonal hematopoiesis of indeterminate potential flagging with curated gene lists

All thresholds are deterministic, frozen from an independent development cohort, and fully auditable.

## Key Features

- **Zero dependencies**: Core classifier uses only the Python standard library (Python 3.7+)
- **Interpretable**: Every tier assignment is traceable to specific evidence layers
- **Graduated confidence**: Four tiers for different downstream applications (clinical, discovery, screening)
- **CHIP-aware**: Integrated clonal hematopoiesis detection prevents false somatic calls from blood-derived variants
- **Caller-agnostic**: Works with Mutect2, VarDict, or any TSV/MAF-producing pipeline
- **Degraded mode**: Operates with reduced layers when PureCN or Mutect2 annotations are unavailable

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

## Usage

### Full mode (Mutect2 + PureCN)

```bash
python slt_classify.py --input annotated_variants.tsv --output slt_output.tsv
```

### Degraded mode (MAF-level, no PureCN)

```bash
python slt_classify.py --input variants.maf --output slt_output.tsv --degraded
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
| `slt_mode` | str | full / degraded |

## Tier Definitions

### Classification Cascade

Conditions are evaluated in order; the first match determines the tier:

| Tier | Meaning | Conditions |
|------|---------|------------|
| **SLT-A** | High confidence somatic | Posterior >= 0.8 AND evidence in {high, medium} AND not chip_likely |
| **SLT-B** | Likely somatic | Posterior >= 0.5, OR (high evidence AND >= 3 layers); chip_likely blocked |
| **SLT-C** | Possible somatic | Posterior >= 0.2, OR evidence in {high, medium} |
| **SLT-D** | Unlikely somatic | All remaining variants |

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

`chip_likely` status blocks both SLT-A and SLT-B assignment, preventing CHIP variants from reaching the most actionable tiers. `chip_likely` variants are downgraded to SLT-C (maximum) or SLT-D.

### Degraded Mode

When PureCN is unavailable (e.g., MAF-level reanalysis without BAM files):

- Layers 1 (POPAF) and 3 (GERMQ) are disabled
- PureCN posterior is set to null
- Maximum achievable tier is SLT-C
- SLT-A and SLT-B are unreachable

## Thresholds

All thresholds are frozen and deterministic:

| Parameter | Value | Source |
|-----------|-------|--------|
| Posterior SLT-A gate | >= 0.8 | PureCN recommended (Riester et al., 2016) |
| Posterior SLT-B gate | >= 0.5 | PureCN recommended |
| Posterior SLT-C gate | >= 0.2 | PureCN recommended |
| POPAF threshold | >= 5.0 | Development cohort |
| gnomAD AF threshold | < 0.001 | Development cohort |
| GERMQ threshold | >= 30 | Development cohort |
| COSMIC confirmed min | >= 5 | Development cohort |
| COSMIC hotspot min | >= 10 | Development cohort |
| CGC + COSMIC sample min | >= 2 | Development cohort |

## Benchmark Performance (SEQC2 HCC1395)

Validated on the SEQC2 HCC1395 breast cancer truth set (455 variants in evaluation regions), processed in tumor-only mode:

| Threshold | Sensitivity | PPV | F1 |
|-----------|-------------|-----|-----|
| >= SLT-A | 18.0% | 78.1% | 0.293 |
| >= SLT-B | 22.4% | 63.7% | 0.332 |
| >= SLT-C | 89.5% | 16.8% | 0.283 |

PureCN posterior AUROC: 0.777 (95% CI: 0.735-0.819).

## Tests

```bash
# Install test dependencies
pip install pytest

# Run tests
pytest tests/ -v
```

The test suite covers all evidence layers, CHIP classification, tier assignment cascade, degraded mode, boundary conditions, and integration tests.

## Repository Structure

```
somatic-likelihood-tiering/
├── slt_classify.py          # Core classifier (standalone, no dependencies)
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

> [Author Names]. Somatic Likelihood Tiering: An Interpretable Multi-Evidence Framework for Somatic Variant Classification in Tumor-Only Whole-Exome Sequencing. *[Journal]*, 2026.

## License

MIT License. See [LICENSE](LICENSE) for details.
