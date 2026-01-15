# A resource for malaria gold-standard validation data

## Aim of this resource
The aim of this resource is to provide datasets that may be used in validating malaria genomics pipelines. Validating a custom pipeline is a key step in the accreditation process.  
Two datasets are provided:  
- high-quality real-world data from public resources
- simulated data: example data and recipes 

## Real-world data
The data in this section was obtained from two high-profile public malaria genomics resources:  
1. [MalariaGEN Pf8](https://www.malariagen.net/data_package/open-dataset-plasmodium-falciparum-v80/)
2. [GenRe Mekong](https://www.malariagen.net/resource/29/)

A large number of samples (see notebooks for details) was examined by both projects under a common sample ID. This provides an opportunity to create high-quality datasets where two different analysis methodologies based on different sequencing technologies arrive at the same conclusions.  

### Assumptions
Analysis pipelines differ in scope and methodology. Thus, no single dataset will be applicable to every possible genomics pipeline. For the selection of real-world public datasets, we are assuming that the pipeline to be validated produces at least some of the following:  
- genotype (haplotype) calls at known drug-resistance loci
- high-level drug-resistance phenotype calls 
Furthermore, we assume that the pipeline can work with _P. falciparum_ data.  

For pipelines that are using amplicon sequencing (AmpSeq) data, we assume that the pipeline can work with the SpotMalaria panel. For details about this panel, consult the [SpotMalaria technical manual](https://ngs.sanger.ac.uk/production/malaria/Resource/29/20200705-GenRe-04a-SpotMalaria-0.39.pdf) and [this SpotMalaria supplementary data file](https://www.malariagen.net/wp-content/uploads/2023/11/20200705-GenRe-04b-SpotMalaria-SupplementaryFile1.xlsx), which provides details for every primer used in the panel.  

For pipelines that work with AmpSeq data for a specific primer panel that is not SpotMalaria, please take a look at the section on [simulated data](#simulated-data), which shows how to create simulated runs for cases where real-world data may not exist.

### Approach for creating high-quality datasets
In order to create the most high-quality datasets, we made use of that fact that a large number of samples have been independently analysed by both projects and that the two projects use different sequencing technology and analysis pipelines to create genotype and resistance phenotype calls. 
To make use of this overlap between the projects, the datasets provided are based on full agreement of results between Pf8 and GenRe Mekong. Thus, every sample in the final datasets have results that are confirmed by two independent sequencing and analysis methods.  

### Datasets provided
1. Samples where all inferred drug resistance phenotypes that have been tested in both projects are identical: [Pf8-GenReMekong_concordant_phenotypes.csv](real-world_gold-standard_data/Pf8-GenReMekong/Pf8-GenReMekong_concordant_phenotypes.csv)
2. All samples where all genotypes at loci known to be relevant to drug resistances are identical: [Pf8-GenReMekong_concordant_genotypes.csv)](real-world_gold-standard_data/Pf8-GenReMekong/Pf8-GenReMekong_concordant_genotypes.csv)
3. A subset of samples, each representing one distinct pattern of drug resistance haplotypes: [Pf8-GenReMekong_concordant_genotypes_representative_samples.csv](real-world_gold-standard_data/Pf8-GenReMekong/Pf8-GenReMekong_concordant_genotypes_representative_samples.csv)

Full details on how these datasets were created from public sources is provided in the form of an executable [Jupyter](https://jupyter.org/) notebook [pf8genre.ipynb](real-world_gold-standard_data/Pf8-GenReMekong/scripts/pf8genre.ipynb). To [run the notebook](https://jupyter.org/install#jupyter-notebook), Python and the pandas library are the only dependencies that need to be installed.  Alternatively, consider using a public service to run Jupyter notebooks such as [Google Colab](https://colab.research.google.com/).

For information about data access, please visit https://malariagen.github.io/parasite-data/pf8/Data_access.html

As detailed in the Jupyter notebook, additional files are available in folder [additional_output_files](real-world_gold-standard_data/Pf8-GenReMekong/additional_output_files). This includes versions of the above files with ENA FTP download links as well as ready-made download manifest data that can be used with a downloader tool provided [here](real-world_gold-standard_data/Pf8-GenReMekong/scripts/ENA_data_helper.py) to retrieve the read FASTQ files for the samples in these datasets.    

### How to use the datasets
The three datasets are provided as comma-separated tables. Most data fields are either directly taken from public data or are calculated from the public data fields as detailed in the accompanying [Jupyter notebook]((real-world_gold-standard_data/Pf8-GenReMekong/scripts/pf8genre.ipynb)). All of the changes to the data columns are limited to renaming columns or extracting values from columns, unchanged, into new columns to enable comparisons between Pf8 and GenRe Mekong data.  Data dictionaries that describe the original public data fields are identified and linked to in the Jupyter notebook.

Each dataset has a column called 'sample', which contains the sample ID that is used in both projects, Pf8 and GenRe Mekong. This is the "primary key" of the data and it can be used to obtain raw sequencing data from public archives. Please note that the field was uploaded as "sample title" to ENA.  

The sample ID can also be used to add more of the original metadata to the datasets, if required. For details on how to do this, consult the [data analysis guides for Pf8](https://malariagen.github.io/parasite-data/pf8/Data_access.html). Note that the sample ID column title in the original Pf8 dataset is spelled with an upper case 'S', whereas the same data column is spelled with a lower case 's' in GenRe Mekong. You may have to convert the title accordingly, depending on where additional metadata is coming from.  

To retrieve the FASTQ files from both projects, a custom [ENA data helper module](real-world_gold-standard_data/Pf8-GenReMekong/scripts/ENA_data_helper.py) is provided alongside the Jupyter notebook. The notebook uses this module to search ENA by sample ID and to add search results back into the sample data tables. At the end of the notebook, a section is provided that creates input files for FASTQ download and demonstrates the use of the ENA data helper on the commandline for the purpose of retrieving the FASTQ files along with a manifest file.  

## Simulated data
While real-world data from public resources are an important part of any validation strategy, such data suffer from two main issues in the context of pipeline validation:  
- data is generated using specific lab techniques that may not be be compatible with a given pipeline
- inferred data such as drug-resistance phenotype, is generated by other pipelines. Even if those pipelines are published and highly trusted, the ground truth cannot be verified

To address these issues, we are providing tools and recipes for the creation of simulated data, as well as some example datasets.  
Simulated data can be made to match the expected input of a given pipeline. For example, a pipeline that expects AmpSeq data that was created using a specific amplicon primer panel should be validated against sequencing data obtained by using that same primer panel. High-quality reference data may not be available in such a case. With simulation, the right type of data can be constructed and known genotypes can be constructed to verify the results of the pipeline run against a ground truth dataset.

### Recipe and tools for creating simulated dataset
Elsewhere in this repository, we provide a [read simulation pipeline](https://github.com/advISO-project/pop_var_sim) that builds on published tools to facilitate the creation of simulated read datasets with known genotypes and the ability to simulate custom AmpSeq panels.  

A [Jupyter](https://jupyter.org/) notebook is provided [here](simulated_data/prepare_simulation_run/scripts/prepare_files_for_simulation.ipynb). It contains a working recipe for creating a simulated read dataset. The recipe starts by using real-world data to obtain realistic guides for read counts to be simulated and provides detailed instructions on how to produce the configuration files that are required by the [read simulation pipeline](https://github.com/advISO-project/pop_var_sim) from information that is publicly available.  

### How to use the simulated data
The notebook includes a command that the user can run locally to produce simulated read FASTQ files with the configuration provided by the recipe in the notebook. This creates a specific set of FASTQ files for WT and a mutant genotype, based on the read counts obtained from high-quality real-world data. While this may be useful as-is, the main aim of the simulation recipe is to provide a template for users to build their own recipes and datasets, tailor-made for their specific pipeline needs.   

For example, the recipe shows how to use a published set of amplicon primer positions to simulate a run with that primer panel. In order to transfer this to a different panel, primer positions need to be obtained from different sources. This may require mapping of primer sequences to the reference genome, depending on the data available for the specific panel. The recipe also shows how a mutant that is provided in the form of an amino acid position in a named gene can be transformed into the vcf-like format of genomic positions and SNP details required by the simulation pipeline.  