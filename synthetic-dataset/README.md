# Synthetic SpotMalaria dataset construction

This directory contains the design and construction materials for the synthetic
SpotMalaria malaria validation dataset. The intended publication dataset contains
25 designed samples with resistance-associated amino-acid states and resistance
classifications expected under the documented interpretation rules. The
classifications are not experimentally observed phenotypes.

The principal construction account is the executable
[`prepare_simulation.ipynb`](notebook/prepare_simulation.ipynb) notebook. The
[`pop_var_sim_input_files`](pop_var_sim_input_files/) directory contains the
simulation inputs, while
[`malaria_DR.synthetic_samples_design.csv`](malaria_DR.synthetic_samples_design.csv)
records the designed states and expected classifications. The accompanying
INSDC manifest associates each sample with its GRC1, GRC2, and SPEC SpotMalaria
read runs.

From the repository root, prepare the machine-generated portion of a candidate
publication bundle with:

```text
lib/publish_synthetic_dataset.py --output-dir <output-directory>
```

The publication script converts every amino-acid call to a canonical JSON set.
Mixed calls record the set of states without encoding mixture proportions or
phase. It also normalises the classification vocabulary, adds the three
associated run accessions, and checks the expected source counts and
relationships before writing output.
