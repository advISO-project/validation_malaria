# Pf8/GenRe concordance-based dataset construction

This directory contains the construction outputs for the real-world malaria
validation dataset. The intended publication dataset contains 74 samples: one
deterministically selected sample for every distinct observed pattern across
the included genotype calls. All are drawn from an eligible population for
which Pf8 and GenRe Mekong report concordant genotype-derived resistance
classifications and concordant, non-missing calls at the selected amino-acid
markers.

The principal construction account is the executable
[`pf8genre.ipynb`](notebook/pf8genre.ipynb) notebook. It produces:

- a 904-record intermediate classification-concordant set;
- the 836-record genotype-concordant eligible population; and
- the 74-record publication set containing the
  lexicographically first sample for each distinct pattern of included genotype
  calls.

The 904- and 836-record sets are construction artefacts, not separate
publication datasets. The 74-record publication set is designed for genotype
pattern coverage, not population or technical representativeness. It does not
preserve genotype-pattern frequencies or repeated observations, and selection
does not deliberately span the range or balance of sequencing coverage, read
quality, country, or AmpSeq process in the eligible population. The sample for
each pattern is selected by sorting sample IDs lexicographically and retaining
the first.

From the repository root, prepare the machine-generated portion of a candidate
publication bundle with:

```text
lib/publish_real_world_dataset.py --output-dir <output-directory>
```

The publication script selects the 74 publication records from the 836-record
eligible genotype table, joins classifications from the 904-record intermediate
table, converts genotype calls to canonical JSON sets, and associates every
record with one Pf8 WGS run and three GenRe SpotMalaria runs. It verifies that
the published records cover all distinct included genotype patterns and that
each is the lexicographically first eligible sample for its pattern before
writing output.
