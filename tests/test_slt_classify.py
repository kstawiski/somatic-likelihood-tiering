#!/usr/bin/env python3
"""Unit tests for SLT classifier.

Run: pytest tests/ -v
"""

import math
import pytest
from slt_classify import (
    safe_float,
    safe_int,
    evaluate_layer1_popaf,
    evaluate_layer2_gnomad,
    evaluate_layer3_germq,
    evaluate_layer4_cosmic,
    count_somatic_layers,
    assign_evidence_level,
    classify_chip,
    assign_slt_tier,
    classify_variant,
    POSTERIOR_A, POSTERIOR_B, POSTERIOR_C,
    POPAF_THRESHOLD, GNOMAD_AF_THRESHOLD, GERMQ_THRESHOLD,
    COSMIC_CONFIRMED_MIN, COSMIC_HOTSPOT_MIN, COSMIC_CGC_SAMPLE_MIN,
)


# =============================================================================
# Parsing helpers
# =============================================================================

class TestSafeFloat:
    def test_valid(self):
        assert safe_float("3.14") == 3.14

    def test_zero(self):
        assert safe_float("0") == 0.0

    def test_none(self):
        assert safe_float(None) is None

    def test_empty(self):
        assert safe_float("") is None

    def test_dot(self):
        assert safe_float(".") is None

    def test_na(self):
        assert safe_float("NA") is None

    def test_nan_string(self):
        assert safe_float("nan") is None

    def test_nan_float(self):
        assert safe_float(float("nan")) is None

    def test_default(self):
        assert safe_float("", default=0.0) == 0.0

    def test_invalid(self):
        assert safe_float("abc") is None


class TestSafeInt:
    def test_valid(self):
        assert safe_int("42") == 42

    def test_float_string(self):
        assert safe_int("3.9") == 3  # truncates

    def test_none(self):
        assert safe_int(None) is None

    def test_empty(self):
        assert safe_int("") is None

    def test_default(self):
        assert safe_int("", default=0) == 0


# =============================================================================
# Evidence Layer Tests
# =============================================================================

class TestLayer1Popaf:
    def test_somatic_high_popaf(self):
        """POPAF >= 5.0 → somatic-supporting."""
        assert evaluate_layer1_popaf(5.0) is True
        assert evaluate_layer1_popaf(40.0) is True

    def test_germline_low_popaf(self):
        """POPAF < 5.0 → NOT somatic-supporting."""
        assert evaluate_layer1_popaf(4.9) is False
        assert evaluate_layer1_popaf(0.0) is False

    def test_missing(self):
        """Missing POPAF → NOT somatic-supporting (by design)."""
        assert evaluate_layer1_popaf(None) is False


class TestLayer2Gnomad:
    def test_somatic_low_af(self):
        """gnomAD_AF < 0.001 → somatic-supporting."""
        assert evaluate_layer2_gnomad(0.0) is True
        assert evaluate_layer2_gnomad(0.0005) is True
        assert evaluate_layer2_gnomad(0.0009) is True

    def test_germline_high_af(self):
        """gnomAD_AF >= 0.001 → NOT somatic-supporting."""
        assert evaluate_layer2_gnomad(0.001) is False
        assert evaluate_layer2_gnomad(0.5) is False

    def test_missing_is_somatic(self):
        """Missing gnomAD_AF → somatic-supporting (by design)."""
        assert evaluate_layer2_gnomad(None) is True


class TestLayer3Germq:
    def test_somatic_high_germq(self):
        """GERMQ >= 30 → somatic-supporting."""
        assert evaluate_layer3_germq(30) is True
        assert evaluate_layer3_germq(100) is True

    def test_germline_low_germq(self):
        """GERMQ < 30 → NOT somatic-supporting."""
        assert evaluate_layer3_germq(29) is False
        assert evaluate_layer3_germq(0) is False

    def test_missing(self):
        """Missing GERMQ → NOT somatic-supporting."""
        assert evaluate_layer3_germq(None) is False


class TestLayer4Cosmic:
    def test_confirmed_somatic(self):
        """COSMIC_CONFIRMED_SOMATIC >= 5 → somatic-supporting."""
        assert evaluate_layer4_cosmic(cosmic_confirmed=5) is True
        assert evaluate_layer4_cosmic(cosmic_confirmed=100) is True

    def test_confirmed_below(self):
        assert evaluate_layer4_cosmic(cosmic_confirmed=4) is False

    def test_hotspot(self):
        """COSMIC_HOTSPOT >= 10 → somatic-supporting."""
        assert evaluate_layer4_cosmic(cosmic_hotspot=10) is True

    def test_hotspot_below(self):
        assert evaluate_layer4_cosmic(cosmic_hotspot=9) is False

    def test_cgc_gene(self):
        """CGC gene + COSMIC_SAMPLE >= 2 → somatic-supporting."""
        assert evaluate_layer4_cosmic(is_cgc_gene=True, cosmic_sample=2) is True

    def test_cgc_gene_below(self):
        assert evaluate_layer4_cosmic(is_cgc_gene=True, cosmic_sample=1) is False

    def test_cgc_false(self):
        assert evaluate_layer4_cosmic(is_cgc_gene=False, cosmic_sample=10) is False

    def test_all_missing(self):
        assert evaluate_layer4_cosmic() is False


# =============================================================================
# Evidence Level Tests
# =============================================================================

class TestEvidenceLevel:
    def test_high(self):
        """high: >= 3 layers AND COSMIC active."""
        assert assign_evidence_level(3, True) == "high"
        assert assign_evidence_level(4, True) == "high"

    def test_medium_3_no_cosmic(self):
        """medium: >= 3 layers but no COSMIC."""
        assert assign_evidence_level(3, False) == "medium"

    def test_medium_2(self):
        """medium: 2 layers."""
        assert assign_evidence_level(2, True) == "medium"
        assert assign_evidence_level(2, False) == "medium"

    def test_low(self):
        """low: < 2 layers."""
        assert assign_evidence_level(1, False) == "low"
        assert assign_evidence_level(0, False) == "low"
        assert assign_evidence_level(1, True) == "low"


# =============================================================================
# CHIP Tests
# =============================================================================

class TestChipClassification:
    def test_tier1_hotspot(self):
        """DNMT3A R882H → chip_likely."""
        assert classify_chip("DNMT3A", 0.05, "R882H", "Missense_Mutation") == "chip_likely"

    def test_tier1_low_vaf(self):
        """Tier1 gene + VAF < 0.15 → chip_likely."""
        assert classify_chip("TET2", 0.10, "Q100X", "Nonsense_Mutation") == "chip_likely"

    def test_tier1_lof_truncating(self):
        """TET2 truncating → chip_likely."""
        assert classify_chip("TET2", 0.30, "Q100X", "Nonsense_Mutation") == "chip_likely"

    def test_tier1_no_criteria(self):
        """Tier1 gene but no specific criteria met → chip_possible."""
        assert classify_chip("DNMT3A", 0.30, "A100V", "Missense_Mutation") == "chip_possible"

    def test_tier2_tp53_hotspot_low_vaf(self):
        """TP53 R175H + low VAF → chip_likely."""
        assert classify_chip("TP53", 0.10, "R175H", "Missense_Mutation") == "chip_likely"

    def test_tier2_tp53_hotspot_high_vaf(self):
        """TP53 R175H + high VAF → chip_possible (not likely)."""
        assert classify_chip("TP53", 0.30, "R175H", "Missense_Mutation") == "chip_possible"

    def test_tier2_kras_hotspot(self):
        """KRAS G12D + low VAF → chip_likely."""
        assert classify_chip("KRAS", 0.05, "G12D", "Missense_Mutation") == "chip_likely"

    def test_tier2_kras_high_vaf(self):
        """KRAS G12D + high VAF → chip_possible."""
        assert classify_chip("KRAS", 0.20, "G12D", "Missense_Mutation") == "chip_possible"

    def test_non_chip_gene(self):
        """Non-CHIP gene → no_chip."""
        assert classify_chip("FGFR3", 0.10, "S249C", "Missense_Mutation") == "no_chip"

    def test_empty_gene(self):
        assert classify_chip("", 0.10, "", "") == "no_chip"

    def test_none_vaf(self):
        """Tier1 gene with missing VAF: hotspot still triggers chip_likely."""
        assert classify_chip("JAK2", None, "V617F", "Missense_Mutation") == "chip_likely"


# =============================================================================
# Tier Assignment Tests
# =============================================================================

class TestSltTier:
    def test_slt_a(self):
        """SLT-A: posterior >= 0.8 AND evidence in {high, medium} AND chip != chip_likely."""
        assert assign_slt_tier(0.95, "high", 4, "no_chip") == "SLT-A"
        assert assign_slt_tier(0.80, "medium", 2, "chip_possible") == "SLT-A"

    def test_slt_a_blocked_by_chip(self):
        """SLT-A blocked by chip_likely → falls to SLT-B."""
        assert assign_slt_tier(0.95, "high", 4, "chip_likely") == "SLT-B"

    def test_slt_a_blocked_by_low_evidence(self):
        """SLT-A blocked by low evidence → falls to SLT-B (posterior >= 0.5)."""
        assert assign_slt_tier(0.95, "low", 1, "no_chip") == "SLT-B"

    def test_slt_b_posterior(self):
        """SLT-B: posterior >= 0.5."""
        assert assign_slt_tier(0.50, "low", 0, "no_chip") == "SLT-B"
        assert assign_slt_tier(0.79, "low", 0, "no_chip") == "SLT-B"

    def test_slt_b_evidence(self):
        """SLT-B: high evidence + >= 3 layers (even with low posterior)."""
        assert assign_slt_tier(0.10, "high", 3, "no_chip") == "SLT-B"

    def test_slt_c_posterior(self):
        """SLT-C: posterior >= 0.2."""
        assert assign_slt_tier(0.20, "low", 0, "no_chip") == "SLT-C"
        assert assign_slt_tier(0.49, "low", 0, "no_chip") == "SLT-C"

    def test_slt_c_evidence(self):
        """SLT-C: medium evidence (even with low posterior)."""
        assert assign_slt_tier(0.10, "medium", 2, "no_chip") == "SLT-C"

    def test_slt_d(self):
        """SLT-D: everything else."""
        assert assign_slt_tier(0.10, "low", 0, "no_chip") == "SLT-D"
        assert assign_slt_tier(0.19, "low", 1, "no_chip") == "SLT-D"
        assert assign_slt_tier(None, "low", 0, "no_chip") == "SLT-D"

    def test_missing_posterior(self):
        """Missing posterior → treated as 0.0 (by design)."""
        assert assign_slt_tier(None, "low", 0, "no_chip") == "SLT-D"
        assert assign_slt_tier(None, "high", 3, "no_chip") == "SLT-B"  # evidence rescue
        assert assign_slt_tier(None, "medium", 2, "no_chip") == "SLT-C"


# =============================================================================
# Degraded Mode Tests
# =============================================================================

class TestDegradedMode:
    def test_degraded_disables_layers(self):
        """Degraded mode: layers 1 and 3 disabled, posterior = None."""
        row = {
            "POPAF": "40",
            "gnomAD_AF": "0.0001",
            "GERMQ": "50",
            "COSMIC_CONFIRMED_SOMATIC": "10",
            "Hugo_Symbol": "TP53",
            "Variant_Classification": "Missense_Mutation",
        }
        result = classify_variant(row, degraded=True)
        assert result["slt_layer1_popaf"] is False  # disabled
        assert result["slt_layer3_germq"] is False  # disabled
        assert result["slt_posterior_used"] == "NA"
        assert result["slt_mode"] == "degraded"

    def test_degraded_keeps_gnomad_cosmic(self):
        """Degraded mode: layers 2 and 4 still active."""
        row = {
            "gnomAD_AF": "0.0001",
            "COSMIC_CONFIRMED_SOMATIC": "10",
            "Hugo_Symbol": "TP53",
            "Variant_Classification": "Missense_Mutation",
        }
        result = classify_variant(row, degraded=True)
        assert result["slt_layer2_gnomad"] is True
        assert result["slt_layer4_cosmic"] is True

    def test_degraded_slt_a_unreachable(self):
        """In degraded mode, SLT-A is unreachable (posterior = 0.0)."""
        row = {
            "gnomAD_AF": "0.00001",
            "COSMIC_CONFIRMED_SOMATIC": "50",
            "COSMIC_HOTSPOT": "20",
            "CGC_GENE": "TRUE",
            "COSMIC_SAMPLE": "100",
            "Hugo_Symbol": "TP53",
            "Variant_Classification": "Missense_Mutation",
        }
        result = classify_variant(row, degraded=True)
        # With 2 layers (gnomAD + COSMIC) and posterior=0.0:
        # evidence = medium (2 layers, not enough for high without layer4+3 layers)
        # SLT-A requires posterior >= 0.8 → unreachable
        assert result["slt_tier"] != "SLT-A"

    def test_degraded_slt_b_unreachable(self):
        """In degraded mode, SLT-B via posterior is unreachable.
        SLT-B via evidence requires high + >= 3 layers, but only 2 available."""
        row = {
            "gnomAD_AF": "0.00001",
            "COSMIC_CONFIRMED_SOMATIC": "10",
            "Hugo_Symbol": "FGFR3",
            "Variant_Classification": "Missense_Mutation",
        }
        result = classify_variant(row, degraded=True)
        # 2 active layers (gnomAD + COSMIC), evidence = medium
        # SLT-B needs posterior >= 0.5 (unavailable) or high + >= 3 (only 2 layers)
        assert result["slt_tier"] in ("SLT-C", "SLT-D")


# =============================================================================
# Full Pipeline Integration Tests
# =============================================================================

class TestFullPipeline:
    def test_high_confidence_somatic(self):
        """Variant with strong somatic evidence → SLT-A."""
        row = {
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
        result = classify_variant(row)
        assert result["slt_tier"] == "SLT-A"
        assert result["slt_n_somatic_layers"] == 4
        assert result["slt_evidence_level"] == "high"
        assert result["slt_chip_status"] == "no_chip"

    def test_germline_variant(self):
        """Variant with strong germline evidence → SLT-D."""
        row = {
            "POSTERIOR.SOMATIC": "0.01",
            "POPAF": "0.5",
            "gnomAD_AF": "0.1",
            "GERMQ": "5",
            "COSMIC_CONFIRMED_SOMATIC": "0",
            "Hugo_Symbol": "BRCA2",
            "t_alt_freq": "0.50",
            "Variant_Classification": "Missense_Mutation",
        }
        result = classify_variant(row)
        assert result["slt_tier"] == "SLT-D"

    def test_chip_variant_downgraded(self):
        """DNMT3A R882H with low VAF → chip_likely blocks SLT-A."""
        row = {
            "POSTERIOR.SOMATIC": "0.90",
            "POPAF": "10",
            "gnomAD_AF": "0.00001",
            "GERMQ": "40",
            "COSMIC_CONFIRMED_SOMATIC": "50",
            "Hugo_Symbol": "DNMT3A",
            "t_alt_freq": "0.05",
            "HGVSp_Short": "R882H",
            "Variant_Classification": "Missense_Mutation",
        }
        result = classify_variant(row)
        assert result["slt_chip_status"] == "chip_likely"
        assert result["slt_tier"] == "SLT-B"  # downgraded from A due to CHIP

    def test_missing_all_fields(self):
        """Empty row → SLT-D with no evidence."""
        row = {}
        result = classify_variant(row)
        assert result["slt_tier"] == "SLT-D"
        assert result["slt_evidence_level"] == "low"
        assert result["slt_chip_status"] == "no_chip"

    def test_full_mode_label(self):
        row = {"POSTERIOR.SOMATIC": "0.5"}
        result = classify_variant(row)
        assert result["slt_mode"] == "full"


# =============================================================================
# Threshold boundary tests (exact values at boundaries)
# =============================================================================

class TestBoundaries:
    def test_posterior_exact_08(self):
        """Posterior exactly 0.8 qualifies for SLT-A gate."""
        assert assign_slt_tier(0.8, "medium", 2, "no_chip") == "SLT-A"

    def test_posterior_just_below_08(self):
        """Posterior 0.799... does NOT qualify for SLT-A."""
        assert assign_slt_tier(0.799, "medium", 2, "no_chip") == "SLT-B"

    def test_posterior_exact_05(self):
        """Posterior exactly 0.5 qualifies for SLT-B."""
        assert assign_slt_tier(0.5, "low", 1, "no_chip") == "SLT-B"

    def test_posterior_exact_02(self):
        """Posterior exactly 0.2 qualifies for SLT-C."""
        assert assign_slt_tier(0.2, "low", 1, "no_chip") == "SLT-C"

    def test_gnomad_exact_threshold(self):
        """gnomAD_AF exactly 0.001 → NOT somatic-supporting."""
        assert evaluate_layer2_gnomad(0.001) is False

    def test_gnomad_just_below(self):
        """gnomAD_AF 0.000999 → somatic-supporting."""
        assert evaluate_layer2_gnomad(0.000999) is True

    def test_popaf_exact_5(self):
        assert evaluate_layer1_popaf(5.0) is True

    def test_germq_exact_30(self):
        assert evaluate_layer3_germq(30) is True


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
