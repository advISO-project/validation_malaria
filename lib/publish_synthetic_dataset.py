#!/usr/bin/env python3
"""Prepare the synthetic SpotMalaria AdvISO dataset for publication."""

from __future__ import annotations

import argparse
from collections import OrderedDict
from pathlib import Path
from typing import Any, Sequence

from publish_dataset_helpers import (
    build_genomic_records,
    build_manifest,
    index_unique,
    normalise_classification,
    parse_state_set,
    read_csv,
    require_count,
    scalar_field,
    state_set_field,
    write_bundle,
)


DATASET_ID = "org.adviso.malaria.synthetic-spotmalaria-drug-resistance"
EXPECTED_RECORDS = 25

GENOTYPE_COLUMNS = OrderedDict(
    [
        ("dhfr_51", "dhfr_51_amino_acid_states"),
        ("dhfr_59", "dhfr_59_amino_acid_states"),
        ("dhfr_108", "dhfr_108_amino_acid_states"),
        ("dhfr_164", "dhfr_164_amino_acid_states"),
        ("crt_76", "crt_76_amino_acid_states"),
        ("dhps_437", "dhps_437_amino_acid_states"),
        ("dhps_540", "dhps_540_amino_acid_states"),
        ("dhps_581", "dhps_581_amino_acid_states"),
        ("dhps_613", "dhps_613_amino_acid_states"),
        ("kelch13_458", "kelch13_458_amino_acid_states"),
        ("kelch13_493", "kelch13_493_amino_acid_states"),
        ("kelch13_539", "kelch13_539_amino_acid_states"),
        ("kelch13_543", "kelch13_543_amino_acid_states"),
        ("kelch13_561", "kelch13_561_amino_acid_states"),
        ("kelch13_580", "kelch13_580_amino_acid_states"),
        ("kelch13_252", "kelch13_252_amino_acid_states"),
        ("kelch13_578", "kelch13_578_amino_acid_states"),
    ]
)

CLASSIFICATION_COLUMNS = OrderedDict(
    [
        ("Artemisinin", "artemisinin_resistance_classification"),
        ("Chloroquine", "chloroquine_resistance_classification"),
        ("Pyrimethamine", "pyrimethamine_resistance_classification"),
        ("Sulfadoxine", "sulfadoxine_resistance_classification"),
        ("S-P", "sulfadoxine_pyrimethamine_resistance_classification"),
    ]
)

CLASSIFICATIONS = ["sensitive", "resistant", "undetermined"]


def _field_definitions() -> dict[str, dict[str, Any]]:
    fields: dict[str, dict[str, Any]] = {
        "record_id": scalar_field(
            "record_key",
            "Stable identifier assigned in the synthetic design.",
            "string",
        )
    }
    for source_name, field_name in CLASSIFICATION_COLUMNS.items():
        fields[field_name] = scalar_field(
            "reference_value",
            f"Expected {source_name} resistance classification under the "
            "documented interpretation rules.",
            "string",
            allowed_values=CLASSIFICATIONS,
        )
    for _, field_name in GENOTYPE_COLUMNS.items():
        marker = field_name.removesuffix("_amino_acid_states")
        gene, residue = marker.rsplit("_", 1)
        fields[field_name] = state_set_field(
            f"Designed {gene} amino-acid states at residue {residue}; "
            "order has no meaning and mixture proportions and phase are not encoded."
        )
    return fields


def build_dataset(
    source_dir: Path,
    *,
    version: str,
    licence: str,
    construction_repository: str,
    construction_revision: str,
    construction_entry_point: str | None,
    contact: str,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Transform the synthetic design and accessions into publication data."""
    design_path = source_dir / "malaria_DR.synthetic_samples_design.csv"
    insdc_path = (
        source_dir
        / "malaria_DR.synthetic_samples_design.INSDC_manifest.spotmalaria.csv"
    )
    required = {"sample", *GENOTYPE_COLUMNS, *CLASSIFICATION_COLUMNS}
    rows = read_csv(design_path, required)
    require_count(len(rows), EXPECTED_RECORDS, "synthetic design records")
    indexed = index_unique(rows, "sample", design_path)

    comparison_records = []
    for sample in sorted(indexed):
        source = indexed[sample]
        record: dict[str, Any] = {"record_id": sample}
        for source_name, field_name in CLASSIFICATION_COLUMNS.items():
            record[field_name] = normalise_classification(source[source_name])
        for source_name, field_name in GENOTYPE_COLUMNS.items():
            record[field_name] = parse_state_set(source[source_name], "/")
        comparison_records.append(record)

    genomic_records = build_genomic_records(
        list(indexed),
        [insdc_path],
        {
            "GRC1": "spotmalaria-grc1",
            "GRC2": "spotmalaria-grc2",
            "SPEC": "spotmalaria-spec",
        },
    )
    manifest = build_manifest(
        dataset_id=DATASET_ID,
        version=version,
        summary=(
            "Synthetic SpotMalaria read sets with designed resistance-associated "
            "genotypes and rule-based expected resistance classifications."
        ),
        licence=licence,
        construction_repository=construction_repository,
        construction_revision=construction_revision,
        construction_path="synthetic-dataset",
        construction_entry_point=construction_entry_point,
        fields=_field_definitions(),
        genomic_records=genomic_records,
        citations=[
            {
                "id": "pf8",
                "citation": (
                    "MalariaGEN et al. Pf8: an open dataset of Plasmodium "
                    "falciparum genome variation in 33,325 worldwide samples. "
                    "Wellcome Open Research. 2025;10:325."
                ),
                "doi": "10.12688/wellcomeopenres.24031.1",
            },
            {
                "id": "who-antimalarial-resistance",
                "citation": (
                    "World Health Organization. Report on antimalarial drug "
                    "efficacy, resistance and response: 10 years of surveillance "
                    "(2010-2019). Geneva: World Health Organization; 2020. "
                    "ISBN 978-92-4-001281-3."
                ),
            },
            {
                "id": "spotmalaria",
                "citation": (
                    "SpotMalaria Platform Technical Notes and Methods, "
                    "version 0.39."
                ),
            },
        ],
        contact=contact,
    )
    return manifest, comparison_records


def _parser() -> argparse.ArgumentParser:
    repository_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-dir",
        type=Path,
        default=repository_root / "synthetic-dataset",
        help="Folder containing the synthetic design and INSDC manifest.",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--version", default="1.0.0")
    parser.add_argument("--licence", default="CC-BY-4.0")
    parser.add_argument(
        "--construction-repository",
        default="https://github.com/advISO-project/validation_malaria",
    )
    parser.add_argument(
        "--construction-revision",
        default=None,
        help="Immutable construction-repository tag or commit.",
    )
    parser.add_argument(
        "--construction-entry-point",
        default="synthetic-dataset/README.md",
        help=(
            "Repository-relative construction overview path; use an empty "
            "value to omit."
        ),
    )
    parser.add_argument(
        "--contact",
        default="https://github.com/advISO-project/validation_datasets/issues",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Replace generated output files."
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the synthetic dataset builder."""
    args = _parser().parse_args(argv)
    manifest, records = build_dataset(
        args.source_dir,
        version=args.version,
        licence=args.licence,
        construction_repository=args.construction_repository,
        construction_revision=(
            args.construction_revision or f"{DATASET_ID}/v{args.version}"
        ),
        construction_entry_point=args.construction_entry_point or None,
        contact=args.contact,
    )
    write_bundle(
        args.output_dir,
        manifest,
        records,
        template_title="Synthetic SpotMalaria drug-resistance comparison",
        template_description=(
            "Template for comparing pipeline results with the canonical "
            "synthetic SpotMalaria reference."
        ),
        overwrite=args.overwrite,
    )
    print(f"Wrote {len(records)} records to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
