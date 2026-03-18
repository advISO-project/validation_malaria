# Methodology and datasets for validating pipelines that predict the drug-resistance status of Plasmodium falciparum samples using sequencing

## Introduction
The aim of this resource is to provide datasets that may be used in validating bioinformatics pipelines that use genomic data to infer the drug resistance status of Plasmodium falciparum samples. Validation of pipelines is key activity in ISO accreditation.

Ideally, the basis for such a dataset would be a set of samples that (a) have been sequenced, and (b) have a drug resistance status (i.e. sensitive or resistent) that has been confirmed in a clinical setting. However, as no such datasets are freely available for malaria, we take two different approaches to fill this gap: one based on real samples for which the genomic data and prediction of the resistance status have been obtained by two independent means, with concordance being used as marker of quality/confidence; and the other based on custom-designed sythetic read sets, for which the "correct answer" is known by construction.

## Approach 1: Real-world data

The data in this approach was obtained from two highly-used publised malaria genomics resources:  
1. Pf8 ([publication](https://pubmed.ncbi.nlm.nih.gov/41267740), and [dataset](https://www.malariagen.net/data_package/open-dataset-plasmodium-falciparum-v80/)). This project has used Illumina whole-genome sequencing (WGS) and a carefully-designed variant characterisation pipeline to genotype over 35,000 Plasmodium falciparum samples.  
2. GenRe Mekong ([publication](https://pubmed.ncbi.nlm.nih.gov/40394107) and [dataset](https://www.malariagen.net/resource/29/)). This project has used the [SpotMalaria](https://ngs.sanger.ac.uk/production/malaria/Resource/29/20200705-GenRe-04a-SpotMalaria-0.39.pdf) to selectively amplify and (Illumina) sequence regions of the genome containing known markers of drug resistance, for just under 10,000 samples.  

Both project perform genotyping and use the genotype calls to infer/predict drug resistance status for a number of drugs. However, the two projects use different sequencing technologies/assays, different genotyping pipelines/workflows, and different processes for inferring DR status using the genotype calls. A large number of samples have been sequenced and analysed by both projects, under a common sample ID. This provides an opportunity to determine where two different approaches have arrived at the same conclusion. 

A dataset comprising the set of samples that agree between the two dataset is useful in a numbers of ways. It represents a set of samples for which the genotype calls and DR resistance status are high-confidence, and thus the statements made should remain true true for other assays used to sequence the same samples (e.g. using a different amplicon panel, or sequencing technology). Furthermore, the (real) sequencing data associated with the samples acts as subset of Pf8 and GeneRe Mekong that can be used to verify the outputs of the existing pipelines used for the projects - or other pipelines designed for the same types of data. 

### Datasets provided
1. Samples where all inferred drug resistance phenotypes that have been tested in both projects are identical: [Pf8-GenReMekong_concordant_phenotypes.csv](real-world-dataset/Pf8-GenReMekong_concordant_phenotypes.csv)
2. The subset of the samples in (1), where, in addtion, all genotypes at loci known to be relevant to drug resistances are identical: [Pf8-GenReMekong_concordant_genotypes.csv](real-world-dataset/Pf8-GenReMekong_concordant_genotypes.csv)
3. A subset of samples in (2), each representing one distinct pattern of drug resistance haplotypes: [Pf8-GenReMekong_concordant_genotypes_representative_samples.csv](real-world-dataset/Pf8-GenReMekong_concordant_genotypes_representative_samples.csv)

Full details on how these datasets were created from public sources is provided in the form of an executable Jupyter notebook [pf8genre.ipynb](real-world-dataset/notebook/pf8genre.ipynb). 

### How to use the datasets
The three datasets are provided as comma-separated tables. Most data fields are either directly taken from public data or are calculated from the public data fields as detailed in the[notebook]((real-world-dataset/notebook/pf8genre.ipynb)). All of the changes to the data columns are limited to renaming columns or extracting values from columns, unchanged, into new columns to enable comparisons between Pf8 and GenRe Mekong data.  Data dictionaries that describe the original public data fields are identified and linked to in the notebook.

Each dataset has a column called 'sample', which contains the sample ID that is used in both projects, Pf8 and GenRe Mekong. This is the "primary key" of the data and it can be used to obtain raw sequencing data from public archives. Please note that the field was uploaded as "sample title" to ENA.  

The sample ID can also be used to add more of the original metadata to the datasets, if required. For details on how to do this, consult the [data analysis guides for Pf8](https://malariagen.github.io/parasite-data/pf8/Data_access.html). Note that the sample ID column title in the original Pf8 dataset is spelled with an upper case 'S', whereas the same data column is spelled with a lower case 's' in GenRe Mekong. You may have to convert the title accordingly, depending on where additional metadata is coming from.  

To retrieve the FASTQ files from both projects, a custom [ENA data helper module](lib/ENA_data_helper.py) is provided. The notebook uses this module to search ENA by sample ID and to add search results back into the sample data tables. At the end of the notebook, a section is provided that creates input files for FASTQ download and demonstrates the use of the ENA data helper on the commandline for the purpose of retrieving the FASTQ files along with a manifest file.  

## Approach 2: Synthetic data
While real-world data from public resources are an important part of any validation strategy, such data suffer from some issues in the context of pipeline validation, such as:  
- Real data may not exist to cover all the different scenarios that you wish to test in your pipeline 
- The associated data was generated using specific lab techniques and assays, which broadly limits their use to pipelines designed to work with that specific type of data
- When developing a pipeline for a novel/new assay (e.g. that based on enrichment of specific genomic loci), high quality real data may not be yet available for validation. 

We have therefore provides some tools and recipes for the creation of designed synthetic data sets (by simulation). To demonstrate the use of these tools, we have applied them to the generation of synthetic validation data set for the aforementioned SpotMalaria panel. 

### Recipe and tools for creating simulated dataset

Our pipeine, [pop_var_sim](https://github.com/advISO-project/pop_var_sim), builds on published tools to facilitate the creation of simulated read datasets with known genotypes and the ability to simulate custom AmpSeq panels.  A Jupyter [notebook](synthetic_data/prepare_simulation_run/scripts/prepare_simulation.ipynb) describes the design of the synthetic data set, and includes code for generating the required input files and configuration for pop_var_sim. 

### Description of the resulting dataset

The dataset has been designed to capture a variety of haplotypes that confer resistance or sensitivity to a number of the drugs used in malaria treatment. This [spreadsheet](synthetic-dataset/malaria_DR.synthetic_samples_design.v1.csv) records the hapolotype and resistance profile for each synthetic sample. A second [spreadsheet](synthetic-dataset/malaria_DR.synthetic_samples_design.v1.INSDC_manifest.csv) records the ENA sample, experiment and run accessions for the three synthetic libraries associated with the SPOT malaria sub-panels (GRC1, GRC2, SPEC).

The data set has been submitted to ENA under a single BioProject [PRJEB109256](https://www.ebi.ac.uk/ena/browser/view/PRJEB109256)/

## Working with the datasets

### Downloading 

There are various ways to download the fastq run data associated with the above datasets (e.g. from NCBI or ENA). We have provided a convenience script that downloads the data sequentially from ENA. This works perfectly well for the SpotMalaria (amplicon) data, but will take much longer for the Pf8 (WGS) data; for that, you may wish to consider other methods if you are impatient. 

The command creates a folder (--out) into which the FASTQ file pairs are written. It also creates a manifest file that lists run accessions, remote FTP URLs and local (absolute) file paths to the downloaded FASTQ files. Example usage for downloading the real-world representative genotypes SpotMalaria dataset: 

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