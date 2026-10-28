#!/usr/bin/env python3
"""Prepare the Pf8/GenRe AdvISO validation dataset for publication."""

from __future__ import annotations

import argparse
from collections import OrderedDict
from pathlib import Path
from typing import Any, Sequence

from publish_dataset_helpers import (
    ConversionError,
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


DATASET_ID = "org.adviso.malaria.pf8-genre-concordant-drug-resistance"
EXPECTED_ELIGIBLE_RECORDS = 836
EXPECTED_INTERMEDIATE_RECORDS = 904
EXPECTED_PUBLICATION_RECORDS = 74

GENOTYPE_COLUMNS = OrderedDict(
    [
        ("core_mutation_dhfr_51", "dhfr_51"),
        ("core_mutation_dhfr_59", "dhfr_59"),
        ("core_mutation_dhfr_108", "dhfr_108"),
        ("core_mutation_dhfr_164", "dhfr_164"),
        ("core_mutation_crt_72", "crt_72"),
        ("core_mutation_crt_74", "crt_74"),
        ("core_mutation_crt_75", "crt_75"),
        ("core_mutation_crt_76", "crt_76"),
        ("core_mutation_crt_326", "crt_326"),
        ("core_mutation_crt_356", "crt_356"),
        ("core_mutation_dhps_436", "dhps_436"),
        ("core_mutation_dhps_437", "dhps_437"),
        ("core_mutation_dhps_540", "dhps_540"),
        ("core_mutation_dhps_581", "dhps_581"),
        ("core_mutation_dhps_613", "dhps_613"),
        ("core_mutation_mdr1_86", "mdr1_86"),
        ("core_mutation_mdr1_184", "mdr1_184"),
        ("core_mutation_mdr1_1246", "mdr1_1246"),
        ("core_mutation_mdr2_484", "mdr2_484"),
        ("core_mutation_fd_193", "fd_193"),
        ("core_mutation_arps10_127", "arps10_127"),
        (
            "core_mutation_kelch13",
            "kelch13_349_726_changes",
        ),
    ]
)

CLASSIFICATION_COLUMNS = OrderedDict(
    [
        ("Artemisinin", "artemisinin"),
        ("Piperaquine", "piperaquine"),
        ("Mefloquine", "mefloquine"),
        ("Chloroquine", "chloroquine"),
        ("Pyrimethamine", "pyrimethamine"),
        ("Sulfadoxine", "sulfadoxine"),
        ("DHA-PPQ", "dha_ppq"),
        ("AS-MQ", "as_mq"),
        ("S-P", "sp"),
        ("S-P-IPTp", "sp_iptp"),
    ]
)

COUNTRY_CODES = {"LA": "Laos", "VN": "Vietnam"}
CLASSIFICATIONS = ["sensitive", "resistant", "undetermined"]


def _source_paths(source_dir: Path) -> dict[str, Path]:
    prefix = "Pf8-GenReMekong_concordant_genotypes"
    return {
        "genotypes": source_dir / f"{prefix}.csv",
        "phenotypes": source_dir
        / "Pf8-GenReMekong_concordant_phenotypes.csv",
        "representatives": source_dir
        / f"{prefix}_representative_samples.csv",
        "wgs": source_dir
        / f"{prefix}_representative_samples.INSDC_manifest.wgs.csv",
        "spotmalaria": source_dir
        / f"{prefix}_representative_samples.INSDC_manifest.spotmalaria.csv",
    }


def _field_definitions() -> dict[str, dict[str, Any]]:
    fields: dict[str, dict[str, Any]] = {
        "record_id": scalar_field(
            "record_key", "Stable source sample identifier.", "string"
        )
    }
    for source_name, field_name in CLASSIFICATION_COLUMNS.items():
        fields[field_name] = scalar_field(
            "reference_value",
            (
                f"Concordant Pf8 and GenRe genotype-derived {source_name} "
                "resistance classification."
            ),
            "string",
            enum_values=CLASSIFICATIONS,
        )
    for _, field_name in GENOTYPE_COLUMNS.items():
        if field_name == "kelch13_349_726_changes":
            fields[field_name] = state_set_field(
                "Concordant nonsynonymous kelch13 amino-acid changes from "
                "residues 349 through 726; order has no meaning. An empty set "
                "means that no nonsynonymous change was reported in the "
                "specified interval (source value WT)."
            )
            continue
        gene, residue = field_name.rsplit("_", 1)
        fields[field_name] = state_set_field(
            f"Concordant {gene} amino acids at residue {residue}; "
            "order has no meaning."
        )
    fields["country"] = scalar_field(
        "context",
        "Country recorded in the source metadata; not a claim about "
        "representativeness or current parasite origin.",
        "string",
        enum_values=["Laos", "Vietnam"],
    )
    fields["ampseq_process"] = scalar_field(
        "context",
        "GenRe Mekong AmpSeq processing version recorded for the sample.",
        "string",
        enum_values=["AmpSeqV1", "AmpSeqV2"],
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
    """Transform current construction outputs into canonical publication data."""
    paths = _source_paths(source_dir)
    genotype_required = {"sample", "country", "ampseq_process", *GENOTYPE_COLUMNS}
    phenotype_required = {
        "sample",
        "country",
        "ampseq_process",
        *CLASSIFICATION_COLUMNS,
    }
    genotype_rows = read_csv(paths["genotypes"], genotype_required)
    phenotype_rows = read_csv(paths["phenotypes"], phenotype_required)
    representative_rows = read_csv(paths["representatives"], genotype_required)
    require_count(
        len(genotype_rows), EXPECTED_ELIGIBLE_RECORDS, "eligible genotype records"
    )
    require_count(
        len(phenotype_rows),
        EXPECTED_INTERMEDIATE_RECORDS,
        "intermediate classification records",
    )
    require_count(
        len(representative_rows),
        EXPECTED_PUBLICATION_RECORDS,
        "genotype-pattern coverage records",
    )

    genotypes = index_unique(genotype_rows, "sample", paths["genotypes"])
    phenotypes = index_unique(phenotype_rows, "sample", paths["phenotypes"])
    representatives = index_unique(
        representative_rows, "sample", paths["representatives"]
    )
    if not set(genotypes) <= set(phenotypes):
        missing = sorted(set(genotypes) - set(phenotypes))
        raise ConversionError(
            "Eligible genotype records are absent from the intermediate set: "
            f"{missing}"
        )
    if not set(representatives) <= set(genotypes):
        missing = sorted(set(representatives) - set(genotypes))
        raise ConversionError(
            f"Coverage records are absent from the eligible set: {missing}"
        )

    def genotype_pattern(row):
        return tuple(row[column] for column in GENOTYPE_COLUMNS)

    pattern_samples: dict[tuple[str, ...], list[str]] = {}
    for sample, row in genotypes.items():
        pattern_samples.setdefault(genotype_pattern(row), []).append(sample)
    require_count(
        len(pattern_samples), EXPECTED_PUBLICATION_RECORDS, "distinct genotype patterns"
    )
    representative_patterns = {
        genotype_pattern(row): sample for sample, row in representatives.items()
    }
    require_count(
        len(representative_patterns),
        EXPECTED_PUBLICATION_RECORDS,
        "distinct publication genotype patterns",
    )
    if set(representative_patterns) != set(pattern_samples):
        raise ConversionError(
            "Publication rows do not cover every distinct eligible "
            "genotype pattern"
        )
    for pattern, samples in pattern_samples.items():
        expected_sample = sorted(samples)[0]
        if representative_patterns[pattern] != expected_sample:
            raise ConversionError(
                "Representative selection is no longer the lexicographically first "
                f"sample for the pattern containing {expected_sample!r}"
            )

    for sample, representative in representatives.items():
        genotype = genotypes[sample]
        for column in genotype_required:
            if representative[column] != genotype[column]:
                raise ConversionError(
                    f"Publication row {sample!r} disagrees with the eligible "
                    f"row in column {column!r}"
                )

    comparison_records = []
    for sample in sorted(representatives):
        genotype = genotypes[sample]
        phenotype = phenotypes[sample]
        try:
            country = COUNTRY_CODES[genotype["country"]]
        except KeyError as error:
            raise ConversionError(
                f"Unexpected country code {genotype['country']!r} for {sample!r}"
            ) from error
        if phenotype["country"] != country:
            raise ConversionError(
                f"Country disagreement for {sample!r}: {genotype['country']!r} "
                f"versus {phenotype['country']!r}"
            )
        if phenotype["ampseq_process"] != genotype["ampseq_process"]:
            raise ConversionError(f"AmpSeq process disagreement for {sample!r}")

        record: dict[str, Any] = {"record_id": sample}
        for source_name, field_name in CLASSIFICATION_COLUMNS.items():
            record[field_name] = normalise_classification(phenotype[source_name])
        for source_name, field_name in GENOTYPE_COLUMNS.items():
            record[field_name] = parse_state_set(
                genotype[source_name],
                ",",
                wild_type_is_empty=(
                    field_name == "kelch13_349_726_changes"
                ),
            )
        record["country"] = country
        record["ampseq_process"] = genotype["ampseq_process"]
        comparison_records.append(record)

    genomic_records = build_genomic_records(
        list(representatives),
        [paths["wgs"], paths["spotmalaria"]],
        {
            "WGS": "pf8-wgs",
            "GRC1": "genre-grc1",
            "GRC2": "genre-grc2",
            "SPEC": "genre-spec",
        },
    )
    manifest = build_manifest(
        dataset_id=DATASET_ID,
        version=version,
        summary=(
            "One deterministically selected Pf8/GenRe Mekong sample for each "
            "distinct observed pattern across the included concordant genotype "
            "calls, with concordant genotype-derived resistance classifications."
        ),
        licence=licence,
        construction_repository=construction_repository,
        construction_revision=construction_revision,
        construction_path="real-world-dataset",
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
                "id": "genre-mekong",
                "citation": (
                    "Genetic surveillance of Plasmodium falciparum populations "
                    "following treatment policy revisions in the Greater Mekong "
                    "Subregion. Nature Communications. 2025;16:4689."
                ),
                "doi": "10.1038/s41467-025-59946-1",
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
        default=repository_root / "real-world-dataset",
        help="Folder containing the real-world construction CSV files.",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--version", default="1.0.0")
    parser.add_argument(
        "--licence", default="LicenseRef-MalariaGEN-Pf8-GenRe"
    )
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
        default="real-world-dataset/README.md",
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
    """Run the real-world dataset builder."""
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
        template_title="Pf8/GenRe concordant drug-resistance comparison",
        template_description=(
            "Template for comparing pipeline results with the canonical "
            "Pf8/GenRe concordance-based reference."
        ),
        overwrite=args.overwrite,
    )
    print(f"Wrote {len(records)} records to {args.output_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ConversionError as error:
        raise SystemExit(f"error: {error}") from error
