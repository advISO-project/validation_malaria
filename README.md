# Methodology and datasets for validating pipelines that report drug-resistance markers and classifications for Plasmodium falciparum samples using sequencing

## Introduction
The aim of this resource is to provide datasets that may be used in validating bioinformatics pipelines that use genomic data to report drug-resistance markers and genotype-derived resistance classifications for Plasmodium falciparum samples. Validation of pipelines is a key activity in ISO accreditation.

Ideally, a validation dataset would include sequenced samples whose resistance phenotypes had been established independently through appropriate clinical or laboratory evidence. No suitable openly available malaria dataset of that kind was identified. This resource therefore uses two complementary approaches: real samples for which Pf8 and GenRe Mekong produced concordant genomic calls using distinct assays and analysis workflows, and simulated read sets with genotypes defined by their construction. The real-world concordance provides evidence for the selected reference calls but does not establish them as independently adjudicated biological or clinical truth. The resistance classifications assigned to the synthetic samples are expected classifications under the documented interpretation rules rather than experimentally observed phenotypes.

## Approach 1: Real-world data

The data in this approach was obtained from two highly-used publised malaria genomics resources:  
1. Pf8 ([publication](https://pubmed.ncbi.nlm.nih.gov/41267740), and [dataset](https://www.malariagen.net/data_package/open-dataset-plasmodium-falciparum-v80/)). This project has used Illumina whole-genome sequencing (WGS) and a carefully-designed variant characterisation pipeline to genotype over 35,000 Plasmodium falciparum samples.  
2. GenRe Mekong ([publication](https://pubmed.ncbi.nlm.nih.gov/40394107) and [dataset](https://www.malariagen.net/resource/29/)). This project has used the [SpotMalaria](https://ngs.sanger.ac.uk/production/malaria/Resource/29/20200705-GenRe-04a-SpotMalaria-0.39.pdf) amplicon panel to selectively amplify and (Illumina) sequence regions of the genome containing known markers of drug resistance, for just under 10,000 samples.  

Both projects perform genotyping and use the genotype calls to infer resistance classifications for a number of drugs. However, the two projects use different sequencing assays, different genotyping workflows, and different processes for inferring resistance classifications from the genotype calls. A large number of samples have been sequenced and analysed by both projects under a common sample ID. This provides an opportunity to determine where two distinct assay and analysis routes have arrived at the same conclusion.

Agreement between Pf8 and GenRe Mekong identifies samples for which two distinct assay and analysis routes produced the same genotype-derived resistance classifications and, for the selected subset, the same amino-acid calls at the specified marker positions. These concordant calls can be used as reference values when assessing other pipelines. Concordance increases confidence in the calls but does not demonstrate that they are correct, independent of shared assumptions, or applicable to every assay. The associated sequencing data also provide test inputs for pipelines intended to process Pf8 WGS or GenRe Mekong SpotMalaria data.

### Construction outputs
1. An intermediate set of 904 samples where all genotype-derived resistance classifications are identical between the two projects: [Pf8-GenReMekong_concordant_phenotypes.csv](real-world-dataset/Pf8-GenReMekong_concordant_phenotypes.csv)
2. An eligible population of 836 samples: the subset of (1) for which all selected amino-acid and normalised kelch13 calls are also concordant and non-missing: [Pf8-GenReMekong_concordant_genotypes.csv](real-world-dataset/Pf8-GenReMekong_concordant_genotypes.csv)
3. The 74-sample publication set selected deterministically from (2), with one sample for each distinct pattern of the selected genotype calls: [Pf8-GenReMekong_concordant_genotypes_representative_samples.csv](real-world-dataset/Pf8-GenReMekong_concordant_genotypes_representative_samples.csv)

The 74-record table is the publication target; the 904- and 836-record tables remain construction artefacts that document how it was derived. The published set covers every distinct observed pattern across the 22 included genotype fields, but it is not intended to represent the source population or the technical variation within it. Selecting one sample per pattern removes information about pattern frequencies and repeated observations, and does not deliberately preserve the range or balance of sequencing coverage, read quality, country, or AmpSeq process found among the 836 eligible samples.

Full details on how these datasets were created from public sources is provided in the form of an executable Jupyter notebook [pf8genre.ipynb](real-world-dataset/notebook/pf8genre.ipynb). 

### How to use the datasets
The intermediate datasets are provided as comma-separated tables. Their fields are derived from the public Pf8 and GenRe Mekong resources using the process documented in the [construction notebook](real-world-dataset/notebook/pf8genre.ipynb). This process includes sample filtering, field selection and renaming, reconciliation of resistance-classification vocabularies, normalisation of amino-acid and kelch13 mutation representations, concordance assessment, selection of genotype-pattern representatives, and matching records to public INSDC accessions. The notebook links to the source data dictionaries and documents the transformations applied.

Each table has a column called `sample`, which contains the sample ID used in both Pf8 and GenRe Mekong. This is the primary key of the construction data and can be used to obtain raw sequencing data from public archives; the corresponding INSDC field is `sample_title`.

The sample ID can also be used to add more of the original metadata to the datasets, if required. For details on how to do this, consult the [data analysis guides for Pf8](https://malariagen.github.io/parasite-data/pf8/Data_access.html).  

## Approach 2: Synthetic data
While real-world data from public resources are an important part of any validation strategy, such data suffer from some issues in the context of pipeline validation, such as:  
- Real data may not exist to cover all the different scenarios that you wish to test in your pipeline 
- The associated data was generated using specific lab techniques and assays, which broadly limits their use to pipelines designed to work with that specific type of data
- When developing a pipeline for a novel/new assay (e.g. that based on enrichment of specific genomic loci), high quality real data may not be yet available for validation. 

We have therefore provides a tool and a recipe for the creation of designed synthetic data sets (by simulation). By way of demonstration, we have applied the tool to the generation of synthetic validation data set for the aforementioned SpotMalaria panel. 

### Recipe and tools for creating simulated dataset

Our pipeline, [pop_var_sim](https://github.com/advISO-project/pop_var_sim), builds on published tools to facilitate creation of simulated read datasets whose designed genotypes are explicitly encoded in the simulation inputs, with support for custom amplicon panels. A Jupyter [notebook](synthetic-dataset/notebook/prepare_simulation.ipynb) describes the design of the synthetic dataset and includes code for generating the required input files and configuration for `pop_var_sim`.

### Description of the resulting dataset

The dataset contains designed marker combinations selected to exercise resistance-classification scenarios for several antimalarial drugs. This [design table](synthetic-dataset/malaria_DR.synthetic_samples_design.csv) records the intended amino-acid states and the resistance classifications expected under the documented interpretation rules. These classifications are not experimentally observed phenotypes. A second [table](synthetic-dataset/malaria_DR.synthetic_samples_design.INSDC_manifest.spotmalaria.csv) records the ENA sample, experiment, and run accessions for the three synthetic libraries associated with the SpotMalaria sub-panels (GRC1, GRC2, and SPEC).

The data set has been submitted to ENA under a single BioProject [PRJEB109256](https://www.ebi.ac.uk/ena/browser/view/PRJEB109256)/

## Preparing publication bundles

Two scripts convert the construction outputs into the canonical format used by
`validation_datasets`:

```text
lib/publish_real_world_dataset.py --output-dir <real-output-directory>
lib/publish_synthetic_dataset.py --output-dir <synthetic-output-directory>
```

The scripts generate `dataset_manifest.json`, `comparison_reference.json`,
`comparison-configuration.template.yaml`, and `checksums.sha256`. If the output
directory already contains the human-maintained `README.md`, `CHANGELOG.md`, or
`LICENCE.md`, those files are left unchanged and included in the generated
checksums. Existing generated files are not replaced unless `--overwrite` is
specified.

The publication scripts use the source files in this repository by default,
accept an alternative `--source-dir`, and expose release metadata such as
`--version`, `--licence`, and `--construction-revision` as command-line options.
They fail if record counts, source vocabularies, genotype-pattern selection,
artifact associations, or other agreed construction invariants have changed.
This is intentional: a changed source release must be reviewed before a new
validation dataset version is published. For the real-world dataset, the
publication script emits the 74 genotype-pattern coverage records and uses the
836-record table to verify their eligibility and deterministic selection.

## Working with the datasets

### Downloading 

There are various ways to download the fastq run data associated with the above datasets (e.g. from NCBI or ENA). We have provided a convenience script that downloads the data sequentially from ENA. This works perfectly well for the SpotMalaria (amplicon) data, but will take much longer for the Pf8 (WGS) data, for which the files are much larger; for that, you may wish to consider other methods if you are impatient. 

The command creates a folder (--output_folder) into which the FASTQ file pairs are written. It also creates a manifest file that lists run accessions, remote FTP URLs and local (absolute) file paths to the downloaded FASTQ files. Example usage for downloading the real-world genotype-pattern coverage dataset for SpotMalaria:

```
lib/ENA_data_helper.py download --insdc_manifest Pf8-GenReMekong_concordant_genotypes_representative_samples.INSDC_manifest.spotmalaria.csv --output_folder <fastq_output_folder>
```

### Using the datasets with AmpRecon

The dataset can be used with any bioinformatics pipeline that processes data from the SpotMalaria amplicon panel (and indeed to validate new pipelines designed for SpotMalaria). However, the current reference pipeline for SPOT malaria - [AmpRecon](https://github.com/genomic-surveillance/AmpRecon) - has certain quirks which require some manipulation of the dataset before it can be used. Specifically, we need to:

- Collate the pair of gzipped fastq files for each run into a single interleaved fastq file per run;
- Create a manifest in the specific format required by AmpRecon

We have provided a script to simplify this task:

```
lib/prep_fqs_and_manifest_for_amprecon.py \
    --insdc_manifest synthetic-dataset/malaria_DR.synthetic_samples_design.v1.INSDC_manifest.spot_malaria.csv \
    --insdc_fastq_folder <path_to_where_fastqs_were_download> \
    --output_folder <path_to_where_collated_fastqs_and_AmpRecon_manifest will be written>
```
