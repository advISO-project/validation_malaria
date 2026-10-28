"""Regression tests for the two validation-dataset publication scripts."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from publish_dataset_helpers import ConversionError, parse_state_set
from publish_real_world_dataset import GENOTYPE_COLUMNS as REAL_GENOTYPE_COLUMNS
from publish_real_world_dataset import main as publish_real_world
from publish_synthetic_dataset import main as publish_synthetic


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = REPOSITORY_ROOT.parent


def _read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def _generated_files(path: Path) -> list[str]:
    return sorted(item.name for item in path.iterdir() if item.is_file())


def test_real_world_publisher_converts_current_sources(tmp_path: Path) -> None:
    output = tmp_path / "real"
    assert publish_real_world(["--output-dir", str(output)]) == 0

    records = _read_json(output / "comparison_reference.json")
    manifest = _read_json(output / "dataset_manifest.json")
    assert len(records) == manifest["record_count"] == 74
    assert "genotype_pattern_representative" not in manifest["fields"]
    assert all("genotype_pattern_representative" not in record for record in records)
    assert {record["country"] for record in records} == {"Laos", "Vietnam"}
    assert {record["ampseq_process"] for record in records} == {
        "AmpSeqV1",
        "AmpSeqV2",
    }

    by_id = {record["record_id"]: record for record in records}
    assert by_id["RCN12025"]["artemisinin"] == "resistant"
    assert by_id["RCN12025"]["kelch13_349_726_changes"] == [
        "C580Y"
    ]
    assert any(
        record["kelch13_349_726_changes"] == [] for record in records
    )
    assert any(
        len(record["crt_76"]) == 2 for record in records
    )
    genotype_fields = list(REAL_GENOTYPE_COLUMNS.values())
    assert len(
        {tuple(str(record[name]) for name in genotype_fields) for record in records}
    ) == 74

    genomic = {
        record["record_id"]: record["retrieval"]
        for record in manifest["genomic_data"]["records"]
    }
    assert [item["id"] for item in genomic["RCN12025"]] == [
        "genre-grc1",
        "genre-grc2",
        "genre-spec",
        "pf8-wgs",
    ]
    assert _generated_files(output) == [
        "checksums.sha256",
        "comparison-configuration.template.yaml",
        "comparison_reference.json",
        "dataset_manifest.json",
    ]


def test_synthetic_publisher_converts_mixtures_and_accessions(
    tmp_path: Path,
) -> None:
    output = tmp_path / "synthetic"
    assert publish_synthetic(["--output-dir", str(output)]) == 0

    records = _read_json(output / "comparison_reference.json")
    manifest = _read_json(output / "dataset_manifest.json")
    assert len(records) == manifest["record_count"] == 25
    by_id = {record["record_id"]: record for record in records}
    assert by_id["Pf_Ch_crt_76_het"]["crt_76"] == ["K", "T"]
    assert (
        by_id["Pf_Ch_crt_76_het"]["chloroquine"] == "undetermined"
    )

    genomic = {
        record["record_id"]: record["retrieval"]
        for record in manifest["genomic_data"]["records"]
    }
    assert [item["id"] for item in genomic["Pf_All_mut"]] == [
        "spotmalaria-grc1",
        "spotmalaria-grc2",
        "spotmalaria-spec",
    ]


@pytest.mark.parametrize("publisher", [publish_real_world, publish_synthetic])
def test_generated_manifest_and_completed_template_validate(
    publisher, tmp_path: Path
) -> None:
    jsonschema = pytest.importorskip("jsonschema")
    output = tmp_path / "output"
    publisher(["--output-dir", str(output)])

    manifest_schema = _read_json(
        PROJECT_ROOT
        / "validation_datasets"
        / "schemas"
        / "dataset-manifest-v1.0.schema.json"
    )
    jsonschema.Draft7Validator(
        manifest_schema, format_checker=jsonschema.FormatChecker()
    ).validate(_read_json(output / "dataset_manifest.json"))

    template = yaml.safe_load(
        (output / "comparison-configuration.template.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert template["metadata"]["id"] is None
    assert template["metadata"]["revision"] is None
    for definition in template["fields"]["definitions"].values():
        if definition["source_fields"]["query"] is None:
            definition["source_fields"]["query"] = definition["source_fields"][
                "reference"
            ]
    template["metadata"]["id"] = "org.example.completed"
    template["metadata"]["revision"] = "1.0.0"
    configuration_schema = _read_json(
        PROJECT_ROOT
        / "validation_toolkit"
        / "src"
        / "adviso_validation"
        / "schemas"
        / "comparison-configuration-v1.0.schema.json"
    )
    jsonschema.Draft7Validator(configuration_schema).validate(template)


def test_synthetic_output_is_deterministic(tmp_path: Path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    publish_synthetic(["--output-dir", str(first)])
    publish_synthetic(["--output-dir", str(second)])
    for name in _generated_files(first):
        assert (first / name).read_bytes() == (second / name).read_bytes()


def test_publisher_refuses_to_overwrite_by_default(tmp_path: Path) -> None:
    output = tmp_path / "synthetic"
    publish_synthetic(["--output-dir", str(output)])
    with pytest.raises(ConversionError, match="Refusing to overwrite"):
        publish_synthetic(["--output-dir", str(output)])
    assert (
        publish_synthetic(["--output-dir", str(output), "--overwrite"]) == 0
    )


def test_state_set_conversion_rejects_duplicate_states() -> None:
    with pytest.raises(ConversionError, match="Duplicate state"):
        parse_state_set("N/N", "/")
