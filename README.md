# A resource for malaria gold-standard validation data

## Introduction
The aim of this resource is to provide datasets that may be used in validating bioinformatics pipelines that use genomic data to infer the drug resistance status of Plasmodium falciparum samples. Validation and verification of pipelines are key activities in the ISO accreditation.   

Ideally, the basis for such a dataset would be a set of samples that (a) have been sequenced, and (b) have a drug resistance status (i.e. sensitive or resistent) that has been confirmed in a clinical setting. However, as no such datasets are freely available for malaria, we take two different approaches to fill this gap: one based on real samples for which the genomic data and prediction of the resistance status have been obtained by two independent means, with concordance being used as marker of quality/confidence; and the other based on custom-designed sythetic read sets (for which the "correct answer" is known by construction).

## Approach 1: Real-world data
The data in this approach was obtained from two high-profile public malaria genomics resources:  
1. [MalariaGEN Pf8](https://www.malariagen.net/data_package/open-dataset-plasmodium-falciparum-v80/)
2. [GenRe Mekong](https://www.malariagen.net/resource/29/)

A large number of samples have been sequenced and analysed by both projects,  under a common sample ID. This provides an opportunity to create high-quality datasets where two different analysis methodologies based on different sequencing technologies arrive at the same conclusions.  

### Assumptions
Analysis pipelines differ in scope and methodology. Thus, no single dataset will be applicable to every possible genomics pipeline. For the selection of real-world public datasets, we are assuming that the pipeline to be validated produces at least some of the following:  
- genotype (haplotype) calls at known drug-resistance loci
- high-level drug-resistance phenotype calls 
Furthermore, we assume that the pipeline can work with _P. falciparum_ data.  

For pipelines that are using amplicon sequencing (AmpSeq) data, we assume that the pipeline can work with the SpotMalaria panel. For details about this panel, consult the [SpotMalaria technical manual](https://ngs.sanger.ac.uk/production/malaria/Resource/29/20200705-GenRe-04a-SpotMalaria-0.39.pdf) and [this SpotMalaria supplementary data file](https://www.malariagen.net/wp-content/uploads/2023/11/20200705-GenRe-04b-SpotMalaria-SupplementaryFile1.xlsx), which provides details for primers used in the panel.  

For pipelines that work with AmpSeq data for a specific primer panel that is not SpotMalaria, please take a look at the section on [simulated data](#simulated-data), which shows how to create simulated runs for cases where real-world data may not exist.

### Datasets provided
1. Samples where all inferred drug resistance phenotypes that have been tested in both projects are identical: [Pf8-GenReMekong_concordant_phenotypes.csv](real-world-dataset/Pf8-GenReMekong_concordant_phenotypes.csv)
2. All samples where all genotypes at loci known to be relevant to drug resistances are identical: [Pf8-GenReMekong_concordant_genotypes.csv](real-world-dataset/Pf8-GenReMekong_concordant_genotypes.csv)
3. A subset of samples, each representing one distinct pattern of drug resistance haplotypes: [Pf8-GenReMekong_concordant_genotypes_representative_samples.csv](real-world-dataset/Pf8-GenReMekong_concordant_genotypes_representative_samples.csv)

Full details on how these datasets were created from public sources is provided in the form of an executable [Jupyter](https://jupyter.org/) notebook [pf8genre.ipynb](real-world-dataset/notebook/pf8genre.ipynb). To [run the notebook](https://jupyter.org/install#jupyter-notebook), Python and the pandas library are the only dependencies that need to be installed.  Alternatively, consider using a public service to run Jupyter notebooks such as [Google Colab](https://colab.research.google.com/).

For information about data access, please visit https://malariagen.github.io/parasite-data/pf8/Data_access.html

As detailed in the Jupyter notebook, additional files are available in [notebook folder](real-world-dataset/notebook). This includes versions of the above files with ENA FTP download links as well as ready-made download manifest data that can be used with a downloader tool provided [here](lib/ENA_data_helper.py) to retrieve the read FASTQ files for the samples in these datasets.    

### How to use the datasets
The three datasets are provided as comma-separated tables. Most data fields are either directly taken from public data or are calculated from the public data fields as detailed in the accompanying [Jupyter notebook]((real-world-dataset/notebook/pf8genre.ipynb)). All of the changes to the data columns are limited to renaming columns or extracting values from columns, unchanged, into new columns to enable comparisons between Pf8 and GenRe Mekong data.  Data dictionaries that describe the original public data fields are identified and linked to in the Jupyter notebook.

Each dataset has a column called 'sample', which contains the sample ID that is used in both projects, Pf8 and GenRe Mekong. This is the "primary key" of the data and it can be used to obtain raw sequencing data from public archives. Please note that the field was uploaded as "sample title" to ENA.  

The sample ID can also be used to add more of the original metadata to the datasets, if required. For details on how to do this, consult the [data analysis guides for Pf8](https://malariagen.github.io/parasite-data/pf8/Data_access.html). Note that the sample ID column title in the original Pf8 dataset is spelled with an upper case 'S', whereas the same data column is spelled with a lower case 's' in GenRe Mekong. You may have to convert the title accordingly, depending on where additional metadata is coming from.  

To retrieve the FASTQ files from both projects, a custom [ENA data helper module](lib/ENA_data_helper.py) is provided alongside the Jupyter notebook. The notebook uses this module to search ENA by sample ID and to add search results back into the sample data tables. At the end of the notebook, a section is provided that creates input files for FASTQ download and demonstrates the use of the ENA data helper on the commandline for the purpose of retrieving the FASTQ files along with a manifest file.  

## Approach 2: Synthetic data
While real-world data from public resources are an important part of any validation strategy, such data suffer from some issues in the context of pipeline validation, such as:  
- Real data may not exist to cover all the different scenarios that you wish to test in your pipeline 
- Data is generated using specific lab techniques that may not be be compatible with a given pipeline
- When developing a pipeline for a specific assay (e.g. that based on enrichment of specific genomic loci), high quality real data may not be yet available for validation. 

We have therefore provides some tools and recipes for the creation of designed synthetic data sets (by simulation). To demonstrate the use of these tools, we have applied them to the generation of synthetic validation data set for the SPOTMalaria panel and [associated pipeline](https://github.com/genomic-surveillance/AmpRecon). This data set has been submitted to ENA under BioProject PRJEB109256. 

### Recipe and tools for creating simulated dataset

We have created a pipeline, [pop_var_sim](https://github.com/advISO-project/pop_var_sim) that builds on published tools to facilitate the creation of simulated read datasets with known genotypes and the ability to simulate custom AmpSeq panels.  

A [Jupyter](https://jupyter.org/) notebook is provided [here](synthetic_data/prepare_simulation_run/scripts/prepare_simulation.ipynb). It describes the design of the synthetic data set, and includes code for generating the required input files and configuration for pop_var_sim. 

### Description of the dataset

As detailed in the notebook, the dataset has been designed to capture a variety of haplotypes that confer resistance or sensitivity to a number of the drugs used in malaria treatment. This [spreadsheet](synthetic-dataset/malaria_DR.synthetic_samples_design.v1.csv) records the hapolotype and resistance profile for each synthetic sample. A second [spreadsheet](synthetic-dataset/malaria_DR.synthetic_samples_design.v1.INSDC_manifest.csv) records the ENA sample, experiment and run accessions for the three synthetic libraries associated with the SPOT malaria sub-panels (GRC1, GRC2, SPEC).

### Notes on using the dataset with AmpRecon

The dataset can be used with any bioinformatics pipeline that processes data from the SPOT Malaria amplicon panel (and indeed to validate new pipelines). However, the current reference pipeline for SPOT malaria, [AmpRecon](https://github.com/genomic-surveillance/AmpRecon) has certain quirk which require some manipulation of the dataset before it can be used. Specifically, we need to:

- Collate the pair of gzipped fastq files for each run into a single interleaved fastq file per run;
- Create a manifest in the specific format required by AmpRecon

We have provided a script to simplify this task:

```
lib/prep_fqs_and_manifest_for_amprecon.py \
    --insdc_manifest synthetic-dataset/malaria_DR.synthetic_samples_design.v1.INSDC_manifest.csv \
    --insdc_fastq_folder <path_to_where_fastqs_where_download> \
    --output_folder <path_to_where_collated_fastqs_and_AmpRecon_manifest will be written>
```