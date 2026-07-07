import pytest

import slt_classify as slt


@pytest.mark.parametrize(
    ("gnomad_state", "cosmic_confirmed", "expected_tier"),
    [
        ("rare_callable", "5", "SLT-C"),
        ("rare_callable", "", "SLT-C"),
        ("common", "5", "SLT-C"),
        ("common", "", "SLT-D"),
        ("unevaluable", "5", "SLT-C"),
        ("unevaluable", "", "SLT-C"),
    ],
)
def test_annotation_only_six_cell_matrix(gnomad_state, cosmic_confirmed, expected_tier):
    row = {
        "gnomAD_state": gnomad_state,
        "COSMIC_CONFIRMED_SOMATIC": cosmic_confirmed,
        "Hugo_Symbol": "BRAF",
    }
    if gnomad_state == "common":
        row["gnomAD_AF"] = "0.01"

    result = slt.classify_variant(row, mode="annotation_only")

    assert result["slt_tier"] == expected_tier
    assert result["slt_gnomad_state"] == gnomad_state


def test_rare_cosmic_negative_diverges_between_annotation_only_and_full_mode():
    row = {
        "gnomAD_state": "rare_callable",
        "COSMIC_CONFIRMED_SOMATIC": "",
        "POSTERIOR.SOMATIC": "0",
        "POPAF": "",
        "GERMQ": "",
        "Hugo_Symbol": "BRAF",
    }

    annotation_only = slt.classify_variant(row, mode="annotation_only")
    full = slt.classify_variant(row, mode="full")

    assert annotation_only["slt_tier"] == "SLT-C"
    assert full["slt_tier"] == "SLT-D"
    assert full["slt_evidence_level"] == "low"


def test_missing_gnomad_annotation_is_unevaluable_not_rare_callable():
    result = slt.classify_variant({"Hugo_Symbol": "BRAF"}, mode="annotation_only")

    assert result["slt_gnomad_state"] == "unevaluable"
    assert result["slt_layer2_gnomad"] is False
    assert result["slt_tier"] == "SLT-C"


def test_low_af_without_callability_context_is_unevaluable():
    result = slt.classify_variant(
        {"gnomAD_AF": "0.00001", "POSTERIOR.SOMATIC": "0", "Hugo_Symbol": "BRAF"},
        mode="full",
    )

    assert result["slt_gnomad_state"] == "unevaluable"
    assert result["slt_layer2_gnomad"] is False


def test_record_observed_flag_without_af_or_an_is_unevaluable():
    row = {"gnomAD_record_observed": "true", "POSTERIOR.SOMATIC": "0", "Hugo_Symbol": "BRAF"}

    result = slt.classify_variant(row, mode="full")

    assert result["slt_gnomad_state"] == "unevaluable"
    assert result["slt_layer2_gnomad"] is False


@pytest.mark.parametrize("field", ["gnomAD_callable", "gnomAD_evaluable"])
def test_callable_flag_without_observed_af_is_unevaluable_in_sites_vcf_mode(field):
    row = {field: "true", "POSTERIOR.SOMATIC": "0", "Hugo_Symbol": "BRAF"}

    result = slt.classify_variant(row, mode="full")

    assert result["slt_gnomad_state"] == "unevaluable"
    assert result["slt_gnomad_state_reason"] == "missing_observed_af_for_sites_vcf"
    assert result["slt_layer2_gnomad"] is False


@pytest.mark.parametrize("field", ["gnomAD_callable", "gnomAD_evaluable"])
def test_callable_flag_without_an_does_not_make_low_af_rare_callable(field):
    row = {"gnomAD_AF": "0.00001", field: "true", "Hugo_Symbol": "BRAF"}
    row["POSTERIOR.SOMATIC"] = "0"

    result = slt.classify_variant(row, mode="full")

    assert result["slt_gnomad_state"] == "unevaluable"
    assert result["slt_layer2_gnomad"] is False


@pytest.mark.parametrize("field", ["gnomAD_callable", "gnomAD_evaluable"])
def test_callable_flag_with_low_an_does_not_make_low_af_rare_callable(field):
    row = {
        "gnomAD_AF": "0.00001",
        "gnomAD_AN": "9999",
        field: "true",
        "POSTERIOR.SOMATIC": "0",
        "Hugo_Symbol": "BRAF",
    }

    result = slt.classify_variant(row, mode="full")

    assert result["slt_gnomad_state"] == "unevaluable"
    assert result["slt_layer2_gnomad"] is False


def test_explicit_common_without_observed_common_af_is_unevaluable():
    row = {"gnomAD_state": "common", "Hugo_Symbol": "BRAF"}

    result = slt.classify_variant(row, mode="annotation_only")

    assert result["slt_gnomad_state"] == "unevaluable"
    assert result["slt_gnomad_state_reason"] == "missing_observed_common_af"
    assert result["slt_tier"] == "SLT-C"


@pytest.mark.parametrize(
    "row",
    [
        {"gnomAD_state": "rare_callable", "gnomAD_AF": "0.01", "gnomAD_AN": "20000"},
        {"gnomAD_state": "common", "gnomAD_AF": "0.00001", "gnomAD_AN": "20000"},
    ],
)
def test_conflicting_explicit_gnomad_state_and_af_an_fails_closed(row):
    row = dict(row)
    row["Hugo_Symbol"] = "BRAF"
    row["POSTERIOR.SOMATIC"] = "0"

    result = slt.classify_variant(row, mode="full")

    assert result["slt_gnomad_state"] == "unevaluable"
    assert result["slt_gnomad_state_reason"] == "conflicting_gnomad_state_af_an"
    assert result["slt_layer2_gnomad"] is False


def test_low_af_with_matched_record_but_no_an_is_unevaluable():
    row = {
        "gnomAD_AF": "0.00001",
        "gnomAD_allele_matched": "true",
        "POSTERIOR.SOMATIC": "0",
        "Hugo_Symbol": "BRAF",
    }

    result = slt.classify_variant(row, mode="full")

    assert result["slt_gnomad_state"] == "unevaluable"
    assert result["slt_layer2_gnomad"] is False


def test_low_af_with_adequate_an_is_rare_callable():
    row = {
        "gnomAD_AF": "0.00001",
        "gnomAD_AN": "20000",
        "POSTERIOR.SOMATIC": "0",
        "Hugo_Symbol": "BRAF",
    }

    result = slt.classify_variant(row, mode="full")

    assert result["slt_gnomad_state"] == "rare_callable"
    assert result["slt_layer2_gnomad"] is True


def test_low_af_below_minimum_an_is_unevaluable():
    row = {
        "gnomAD_AF": "0.00001",
        "gnomAD_AN": "9999",
        "POSTERIOR.SOMATIC": "0",
        "Hugo_Symbol": "BRAF",
    }

    result = slt.classify_variant(row, mode="full")

    assert result["slt_gnomad_state"] == "unevaluable"
    assert result["slt_layer2_gnomad"] is False


def test_common_af_is_common_without_needing_callability_flag():
    row = {"gnomAD_AF": "0.01", "Hugo_Symbol": "BRAF"}

    result = slt.classify_variant(row, mode="annotation_only")

    assert result["slt_gnomad_state"] == "common"
    assert result["slt_layer2_gnomad"] is False
    assert result["slt_tier"] == "SLT-D"


def test_prior_output_slt_gnomad_state_is_not_reused_as_input():
    row = {"slt_gnomad_state": "rare_callable", "POSTERIOR.SOMATIC": "0", "Hugo_Symbol": "BRAF"}

    result = slt.classify_variant(row, mode="full")

    assert result["slt_gnomad_state"] == "unevaluable"
    assert result["slt_gnomad_state_reason"] == "missing_or_insufficient_callability"


def test_high_evidence_full_mode_routes_to_slt_b_without_redundant_layer_check():
    tier = slt.assign_slt_tier(posterior=0.0, evidence_level="high", chip_status="no_chip")

    assert tier == "SLT-B"


def test_full_mode_requires_purecn_somatic_probability_column_presence():
    with pytest.raises(ValueError, match="requires PureCN somatic probability column"):
        slt.classify_variant(
            {"gnomAD_AF": "0.00001", "gnomAD_AN": "20000", "Hugo_Symbol": "BRAF"},
            mode="full",
        )


def test_full_mode_blank_posterior_preserves_submitted_zero_posterior_semantics():
    row_high = {
        "Hugo_Symbol": "BRAF",
        "POSTERIOR.SOMATIC": "",
        "POPAF": "5",
        "GERMQ": "",
        "gnomAD_state": "rare_callable",
        "COSMIC_CONFIRMED_SOMATIC": "5",
    }
    result_high = slt.classify_variant(row_high, mode="full")
    assert result_high["slt_evidence_level"] == "high"
    assert result_high["slt_tier"] == "SLT-B"
    assert result_high["slt_posterior_status"] == "absent"
    assert result_high["slt_posterior_used"] == "NA"

    row_medium = {
        "Hugo_Symbol": "BRAF",
        "POSTERIOR.SOMATIC": "",
        "gnomAD_state": "rare_callable",
        "COSMIC_CONFIRMED_SOMATIC": "5",
    }
    result_medium = slt.classify_variant(row_medium, mode="full")
    assert result_medium["slt_evidence_level"] == "medium"
    assert result_medium["slt_tier"] == "SLT-C"

    row_low = {
        "Hugo_Symbol": "BRAF",
        "POSTERIOR.SOMATIC": "",
        "gnomAD_state": "rare_callable",
    }
    result_low = slt.classify_variant(row_low, mode="full")
    assert result_low["slt_tier"] == "SLT-D"


def test_full_mode_blank_posterior_high_evidence_chip_likely_routes_to_slt_c():
    row = {
        "Hugo_Symbol": "DNMT3A",
        "HGVSp_Short": "p.R882H",
        "POSTERIOR.SOMATIC": "",
        "POPAF": "5",
        "gnomAD_state": "rare_callable",
        "COSMIC_CONFIRMED_SOMATIC": "5",
    }

    result = slt.classify_variant(row, mode="full")

    assert result["slt_evidence_level"] == "high"
    assert result["slt_chip_status"] == "chip_likely"
    assert result["slt_tier"] == "SLT-C"


def test_multiallelic_gnomad_record_matches_candidate_alt_allele():
    candidate = ("chr1", "100", "A", "G")
    gnomad_record = ("1", "100", "A", "C,G")

    assert slt.variant_keys_intersect(candidate, gnomad_record) is True


def test_common_prefix_suffix_variant_normalization_matches_complex_snv():
    left_padded = ("chr2", "200", "ATC", "AGC")
    minimal = ("2", "201", "T", "G")

    assert slt.variant_keys_intersect(left_padded, minimal) is True


def test_vcf_anchored_indel_matches_after_chromosome_normalization():
    candidate = ("chr3", "300", "AT", "A")
    gnomad_record = ("3", "300", "AT", "A,T")

    assert slt.variant_keys_intersect(candidate, gnomad_record) is True


def test_float_formatted_position_normalizes_for_variant_matching():
    candidate = ("chr3", "300.0", "AT", "A")
    gnomad_record = ("3", "300", "AT", "A,T")

    assert slt.variant_keys_intersect(candidate, gnomad_record) is True


def test_full_mode_process_file_fails_on_unevaluable_gnomad(tmp_path):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"
    input_path.write_text(
        "Hugo_Symbol\tgnomAD_AF\tPOSTERIOR.SOMATIC\n"
        "BRAF\t0.00001\t0.1\n"
    )

    with pytest.raises(SystemExit) as exc:
        slt.process_file(str(input_path), str(output_path), mode="full")

    assert exc.value.code == 2
    assert not output_path.exists()


def test_full_mode_process_file_fails_on_all_blank_posterior_column(tmp_path):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"
    input_path.write_text(
        "Hugo_Symbol\tgnomAD_AF\tgnomAD_AN\tPOSTERIOR.SOMATIC\n"
        "BRAF\t0.00001\t20000\t\n"
        "KRAS\t0.00001\t20000\t\n"
    )

    with pytest.raises(SystemExit) as exc:
        slt.process_file(
            str(input_path),
            str(output_path),
            mode="full",
            allow_unevaluable_gnomad=True,
        )

    assert exc.value.code == 2
    assert not output_path.exists()


def test_full_mode_process_file_can_explicitly_allow_unevaluable_gnomad(tmp_path):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"
    input_path.write_text(
        "Hugo_Symbol\tgnomAD_AF\tPOSTERIOR.SOMATIC\n"
        "BRAF\t0.00001\t0.1\n"
    )

    slt.process_file(
        str(input_path),
        str(output_path),
        mode="full",
        allow_unevaluable_gnomad=True,
    )

    text = output_path.read_text()
    assert "slt_schema_version" in text
    assert "unevaluable" in text
