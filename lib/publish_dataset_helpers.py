#!/usr/bin/env python3
"""Shared helpers for publishing AdvISO validation dataset bundles."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


class ConversionError(ValueError):
    """Raised when construction input violates an expected invariant."""


SEMVER_PATTERN = re.compile(
    r"^(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)\."
    r"(?:0|[1-9][0-9]*)(?:-[0-9A-Za-z-]+(?:\.[0-9A-Za-z-]+)*)?"
    r"(?:\+[0-9A-Za-z-]+(?:\.[0-9A-Za-z-]+)*)?$"
)


def read_csv(path: Path, required_columns: Iterable[str]) -> list[dict[str, str]]:
    """Read a CSV file and require named columns and non-empty cells."""
    required = tuple(required_columns)
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                raise ConversionError(f"CSV has no header: {path}")
            missing = sorted(set(required) - set(reader.fieldnames))
            if missing:
                raise ConversionError(
                    f"CSV {path} is missing required columns: {', '.join(missing)}"
                )
            rows = list(reader)
    except FileNotFoundError as error:
        raise ConversionError(f"Required source file does not exist: {path}") from error

    for line_number, row in enumerate(rows, start=2):
        for column in required:
            if row[column] is None or not row[column].strip():
                raise ConversionError(
                    f"CSV {path} has an empty {column!r} value at line "
                    f"{line_number}"
                )
    return rows


def index_unique(
    rows: Iterable[Mapping[str, str]], key: str, source: Path
) -> dict[str, Mapping[str, str]]:
    """Index rows by a required unique key."""
    indexed: dict[str, Mapping[str, str]] = {}
    for row in rows:
        value = row[key]
        if value in indexed:
            raise ConversionError(f"Duplicate {key} {value!r} in {source}")
        indexed[value] = row
    return indexed


def require_count(actual: int, expected: int, description: str) -> None:
    """Require an exact count with a construction-oriented error message."""
    if actual != expected:
        raise ConversionError(
            f"Expected {expected} {description}, but found {actual}; "
            "review the source update before publishing a new dataset version"
        )


def require_equal_sets(
    actual: Iterable[str], expected: Iterable[str], description: str
) -> None:
    """Require two identifier collections to contain exactly the same values."""
    actual_set = set(actual)
    expected_set = set(expected)
    if actual_set == expected_set:
        return
    missing = sorted(expected_set - actual_set)
    unexpected = sorted(actual_set - expected_set)
    details = []
    if missing:
        details.append(f"missing={missing[:10]}")
    if unexpected:
        details.append(f"unexpected={unexpected[:10]}")
    raise ConversionError(f"{description} differ: {'; '.join(details)}")


def parse_state_set(
    value: str, delimiter: str, *, wild_type_is_empty: bool = False
) -> list[str]:
    """Convert a delimited source call into a canonical sorted JSON set."""
    text = value.strip()
    if wild_type_is_empty and text == "WT":
        return []
    states = [part.strip() for part in text.split(delimiter)]
    if any(not state for state in states):
        raise ConversionError(f"Invalid empty state in source value {value!r}")
    if len(states) != len(set(states)):
        raise ConversionError(f"Duplicate state in source value {value!r}")
    return sorted(states)


def normalise_classification(value: str) -> str:
    """Convert a source resistance classification to the canonical vocabulary."""
    canonical = value.strip().lower()
    allowed = {"sensitive", "resistant", "undetermined"}
    if canonical not in allowed:
        raise ConversionError(
            f"Unexpected resistance classification {value!r}; expected one of "
            f"{sorted(allowed)}"
        )
    return canonical


def build_genomic_records(
    record_ids: Sequence[str],
    manifest_paths: Sequence[Path],
    retrieval_ids: Mapping[str, str],
) -> list[dict[str, Any]]:
    """Build deterministic per-record retrieval lists from INSDC manifest CSVs."""
    columns = {
        "sample",
        "sample_accession",
        "library_name",
        "experiment_accession",
        "run_accession",
    }
    by_sample: dict[str, dict[str, str]] = {record_id: {} for record_id in record_ids}
    seen_runs: set[str] = set()

    for path in manifest_paths:
        rows = read_csv(path, columns)
        for row in rows:
            sample = row["sample"]
            if sample not in by_sample:
                raise ConversionError(
                    f"INSDC manifest {path} contains unpublished sample {sample!r}"
                )
            prefix = f"{sample}_"
            library = row["library_name"]
            if not library.startswith(prefix):
                raise ConversionError(
                    f"Library {library!r} does not begin with expected {prefix!r}"
                )
            suffix = library[len(prefix) :]
            try:
                retrieval_id = retrieval_ids[suffix]
            except KeyError as error:
                raise ConversionError(
                    f"Unrecognised library suffix {suffix!r} in {path}"
                ) from error
            if retrieval_id in by_sample[sample]:
                raise ConversionError(
                    f"Sample {sample!r} has duplicate retrieval ID {retrieval_id!r}"
                )
            run_accession = row["run_accession"]
            if run_accession in seen_runs:
                raise ConversionError(
                    f"Run accession {run_accession!r} occurs more than once"
                )
            seen_runs.add(run_accession)
            by_sample[sample][retrieval_id] = run_accession

    expected_retrieval_ids = set(retrieval_ids.values())
    records = []
    for record_id in sorted(record_ids):
        observed = set(by_sample[record_id])
        if observed != expected_retrieval_ids:
            raise ConversionError(
                f"Sample {record_id!r} has retrieval IDs {sorted(observed)}, "
                f"expected {sorted(expected_retrieval_ids)}"
            )
        records.append(
            {
                "record_id": record_id,
                "retrieval": [
                    {
                        "id": retrieval_id,
                        "accession": by_sample[record_id][retrieval_id],
                    }
                    for retrieval_id in sorted(expected_retrieval_ids)
                ],
            }
        )
    return records


def scalar_field(
    role: str,
    description: str,
    logical_type: str,
    *,
    allowed_values: Sequence[Any] | None = None,
) -> dict[str, Any]:
    """Create a required scalar field definition for a dataset manifest."""
    definition: dict[str, Any] = {
        "role": role,
        "description": description,
        "value_type": {"logical_type": logical_type},
        "required": True,
    }
    if allowed_values is not None:
        definition["allowed_values"] = list(allowed_values)
    return definition


def state_set_field(description: str) -> dict[str, Any]:
    """Create a required unordered string-set field definition."""
    return {
        "role": "reference_value",
        "description": description,
        "value_type": {"collection_type": "set", "element_type": "string"},
        "required": True,
    }


def build_manifest(
    *,
    dataset_id: str,
    version: str,
    summary: str,
    licence: str,
    construction_repository: str,
    construction_revision: str,
    construction_path: str,
    construction_entry_point: str | None,
    fields: Mapping[str, Mapping[str, Any]],
    genomic_records: Sequence[Mapping[str, Any]],
    citations: Sequence[Mapping[str, str]],
    contact: str,
) -> dict[str, Any]:
    """Build the common manifest envelope."""
    if not SEMVER_PATTERN.fullmatch(version):
        raise ConversionError(
            f"Dataset version is not semantic versioning: {version!r}"
        )
    construction = {
        "repository": construction_repository,
        "revision": construction_revision,
        "path": construction_path,
    }
    if construction_entry_point:
        construction["entry_point"] = construction_entry_point
    return {
        "dataset_manifest_version": "1.0",
        "dataset": {
            "id": dataset_id,
            "version": version,
            "status": "active",
            "summary": summary,
        },
        "organism": {
            "scientific_name": "Plasmodium falciparum",
            "ncbi_taxonomy_id": 5833,
        },
        "record_count": len(genomic_records),
        "licence": licence,
        "contacts": {"issues": contact},
        "construction": construction,
        "citations": list(citations),
        "fields": dict(fields),
        "genomic_data": {
            "insdc_record_type": "reads",
            "records": list(genomic_records),
        },
    }


def validate_bundle(
    manifest: Mapping[str, Any], comparison_records: Sequence[Mapping[str, Any]]
) -> None:
    """Validate important manifest/reference relationships before writing."""
    if manifest["record_count"] != len(comparison_records):
        raise ConversionError("Manifest record_count does not match comparison records")

    fields = manifest["fields"]
    required_fields = {
        name for name, definition in fields.items() if definition["required"]
    }
    record_ids: list[str] = []
    for record in comparison_records:
        record_id = record.get("record_id")
        if not isinstance(record_id, str) or not record_id:
            raise ConversionError("Every comparison record needs a non-empty record_id")
        record_ids.append(record_id)
        unknown = set(record) - set(fields)
        missing = required_fields - set(record)
        if unknown or missing:
            raise ConversionError(
                f"Record {record_id!r} has field mismatch: "
                f"unknown={sorted(unknown)}, missing={sorted(missing)}"
            )
        for name, value in record.items():
            definition = fields[name]
            value_type = definition["value_type"]
            if "collection_type" in value_type:
                if not isinstance(value, list) or any(
                    not isinstance(item, str) for item in value
                ):
                    raise ConversionError(
                        f"Record {record_id!r} field {name!r} is not a string array"
                    )
                if len(value) != len(set(value)):
                    raise ConversionError(
                        f"Record {record_id!r} field {name!r} has duplicate set values"
                    )
            else:
                logical_type = value_type["logical_type"]
                expected_types = {
                    "string": (str,),
                    "boolean": (bool,),
                    "integer": (int,),
                    "decimal": (int, float),
                }[logical_type]
                if not isinstance(value, expected_types):
                    raise ConversionError(
                        f"Record {record_id!r} field {name!r} has the wrong type"
                    )
                allowed = definition.get("allowed_values")
                if allowed is not None and value not in allowed:
                    raise ConversionError(
                        f"Record {record_id!r} field {name!r} has disallowed value "
                        f"{value!r}"
                    )

    if len(record_ids) != len(set(record_ids)):
        raise ConversionError("Comparison record IDs are not unique")
    genomic_ids = [
        record["record_id"] for record in manifest["genomic_data"]["records"]
    ]
    require_equal_sets(genomic_ids, record_ids, "Genomic and comparison record IDs")
    if len(genomic_ids) != len(set(genomic_ids)):
        raise ConversionError("Genomic record IDs are not unique")


def _yaml_quote(value: str) -> str:
    return json.dumps(value, ensure_ascii=False)


def comparison_template(
    manifest: Mapping[str, Any], title: str, description: str
) -> str:
    """Create a deliberately incomplete comparison-configuration YAML template."""
    fields = manifest["fields"]
    in_scope = [
        name
        for name, definition in fields.items()
        if definition["role"] in {"record_key", "reference_value"}
    ]
    lines = [
        "# AdvISO comparison-configuration template.",
        "# Replace each non-key query null and the metadata nulls before use.",
        "# Add query-side representation, normalisation, and mappings as needed.",
        "# The null placeholders deliberately make this template non-runnable.",
        "",
        'schema_version: "1.0"',
        "",
        "metadata:",
        "  id: null",
        "  revision: null",
        f"  title: {_yaml_quote(title)}",
        f"  description: {_yaml_quote(description)}",
        "",
        "fields:",
        "  record_key:",
        '    - "record_id"',
        "  in_scope:",
    ]
    lines.extend(f"    - {_yaml_quote(name)}" for name in in_scope)
    lines.append("  definitions:")

    for name in in_scope:
        definition = fields[name]
        value_type = definition["value_type"]
        lines.extend(
            [
                f"    {name}:",
                "      source_fields:",
                f"        reference: {_yaml_quote(name)}",
                (
                    f"        query: {_yaml_quote(name)}"
                    if definition["role"] == "record_key"
                    else "        query: null"
                ),
                "      required:",
                "        reference: true",
                "        query: true",
            ]
        )
        if "logical_type" in value_type:
            lines.append(
                f"      logical_type: {_yaml_quote(value_type['logical_type'])}"
            )
            comparator = "typed_exact"
        else:
            lines.extend(
                [
                    "      logical_type:",
                    "        collection_type: "
                    f"{_yaml_quote(value_type['collection_type'])}",
                    f"        element_type: {_yaml_quote(value_type['element_type'])}",
                    "        element_comparator:",
                    '          name: "typed_exact"',
                ]
            )
            comparator = f"{value_type['collection_type']}_exact"
        lines.extend(
            [
                "      null_policy:",
                "        allowed: false",
                "        nulls_equal: false",
                "      comparator:",
                f"        name: {_yaml_quote(comparator)}",
                "",
            ]
        )
    return "\n".join(lines).rstrip() + "\n"


def write_json(path: Path, value: Any) -> None:
    """Write deterministic UTF-8 JSON with a trailing newline."""
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def write_checksums(output_dir: Path) -> None:
    """Write SHA-256 checksums for every regular bundle file except itself."""
    checksum_path = output_dir / "checksums.sha256"
    lines = []
    for path in sorted(output_dir.rglob("*")):
        if not path.is_file() or path == checksum_path:
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        relative = path.relative_to(output_dir).as_posix()
        lines.append(f"{digest}  {relative}")
    checksum_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_bundle(
    output_dir: Path,
    manifest: Mapping[str, Any],
    comparison_records: Sequence[Mapping[str, Any]],
    *,
    template_title: str,
    template_description: str,
    overwrite: bool,
) -> None:
    """Validate and write the generated portion of a dataset bundle."""
    validate_bundle(manifest, comparison_records)
    output_dir.mkdir(parents=True, exist_ok=True)
    generated = [
        output_dir / "dataset_manifest.json",
        output_dir / "comparison_reference.json",
        output_dir / "comparison-configuration.template.yaml",
        output_dir / "checksums.sha256",
    ]
    existing = [path for path in generated if path.exists()]
    if existing and not overwrite:
        raise ConversionError(
            "Refusing to overwrite generated files: "
            + ", ".join(str(path) for path in existing)
        )
    write_json(output_dir / "dataset_manifest.json", manifest)
    write_json(output_dir / "comparison_reference.json", comparison_records)
    (output_dir / "comparison-configuration.template.yaml").write_text(
        comparison_template(manifest, template_title, template_description),
        encoding="utf-8",
    )
    write_checksums(output_dir)
