# A resource for malaria gold-standard validation data

## Aim of this resource
The aim of this resource is to provide datasets that may be used in validating malaria genomics pipelines in order to provide evidence that can be used during the accreditation  process.  
Two datasets are provided:  
- high-quality real-world data from public resources
- simulated data: example data and recipes


## Details about the data sets
This section provides details about the provided datasets, including recipes for creating the data, which may serve as a useful starting point for creating tailored datasets for different purposes.  

### Real-world data
The data in this section was obtained from two high-profile public malaria genomics resources:  
1. [MalariaGEN Pf8](https://www.malariagen.net/data_package/open-dataset-plasmodium-falciparum-v80/)
2. [GenRe Mekong](https://www.malariagen.net/resource/29/)

#### Assumptions
Analysis pipelines differ in scope and methodology. Thus, no single dataset will be applicable to every possible genomics pipeline. For the selection of real-world public datasets, we are assuming that the pipeline to be validated produces:  
- amino acid calls for variant loci that are relevant to drug resistance
and/or  
- high-level drug-resistance phenotype calls 
Furthermore, we assume that the pipeline can work with _P. falciparum_ data.

For pipelines that are using amplicon sequencing data, we assume that the pipeline can work with the SpotMalaria panel. For details about this panel, consult the [SpotMalaria technical manual](https://ngs.sanger.ac.uk/production/malaria/Resource/29/20200705-GenRe-04a-SpotMalaria-0.39.pdf) and [this SpotMalaria supplementary data file](https://www.malariagen.net/wp-content/uploads/2023/11/20200705-GenRe-04b-SpotMalaria-SupplementaryFile1.xlsx), which provides details for every primer used in the panel.

#### Approach for creating high-quality datasets
In order to create the most high-quality datasets, we made use of that fact that some samples have been independently analysed by both projects and that the two projects use different sequencing technology and analysis pipelines to create genotype and resistance phenotype calls. 
To make use of this overlap between the projects, the datasets provided are based on full agreement of results between Pf8 and GenRe Mekong. 

#### Datasets provided
1. Samples where all inferred drug resistance phenotypes that have been tested in both projects are identical: [Pf8-GenReMekong_concordant_phenotype_high_quality_samples.csv](real-world_gold-standard_data/Pf8-GenReMekong/Pf8-GenReMekong_concordant_phenotype_high_quality_samples.csv)
2. All samples where all genotypes at loci known to be relevant to drug resistances are identical: [Pf8-GenReMekong_fully_concordant_genotype_high_quality_samples.csv](real-world_gold-standard_data/Pf8-GenReMekong/Pf8-GenReMekong_fully_concordant_genotype_high_quality_samples.csv)
3. A subset of samples, each representing one distinct pattern of drug resistance haplotypes: [Pf8-GenReMekong_fully_concordant_genotype_patterns_representative_samples.csv](real-world_gold-standard_data/Pf8-GenReMekong/Pf8-GenReMekong_fully_concordant_genotype_patterns_representative_samples.csv)

Full details on how these datasets were created from public sources is provided in the form of a [Jupyter](https://jupyter.org/) notebook [pf8genre.ipynb](real-world_gold-standard_data/Pf8-GenReMekong/pf8genre.ipynb). To [run the notebook](https://jupyter.org/install#jupyter-notebook), Python and the pandas library are the only dependencies that need to be installed.  Alternatively, consider using a public service to run Jupyter notebooks such as [Google Colab](https://colab.research.google.com/).

For information about data access, please visit https://malariagen.github.io/parasite-data/pf8/Data_access.html

#### How to use the datasets
The three datasets are provided as comma-separated tables in [this](real-world_gold-standard_data/Pf8-GenReMekong/) folder. Most data fields are either directly taken from public data or are calculaderived from the public data fields as detailed in the accompanying [Jupyter notebook](real-world_gold-standard_data/Pf8-GenReMekong/pf8genre.ipynb). All of the changes to the data columns are limited to renaming columns or extracting values from columns, unchanged, into new columns to enable comparisons between Pf8 and GenRe Mekong data.  Data dictionaries that describe the original public data fields are identified and linked to in the Jupyter notebook.

Each dataset has a column called 'sample', which contains the sample ID that is used in both projects, Pf8 and GenRe Mekong. This is the "primary key" of the data and it can be used to obtain raw sequencing data from public archives. 

The sample ID can also be used to add more of the original metadata to the datasets, if required. For details on how to do this, consult the [data analysis guides for Pf8](https://malariagen.github.io/parasite-data/pf8/Data_access.html). Note that the sample ID column in the original Pf8 dataset is spelled with an upper case 'S', whereas the same data column is spelled with a lower case 's' in GenRe Mekong. 


### Simulated data
While real-world data from public resources are an important part of any validation strategy, such data suffer from two main issues:  
- data is generated using specific lab techniques that may not be be compatible with a given pipeline
- 
- 
- inferred drug-resistance phenotype, are generated by other pipelines. Even if those pipelines are published and highly trusted, the ground truth cannot be verified

To address these issues, we are providing tools and recipes for the creation of simulated data, as well as some example datasets.  
Simulated data can be made to match the expected input of a given pipeline. For example, a pipeline that expects amplicon sequencing data that was created using a specific amplicon primer panel cannot be validated against whole-genome sequencing data. But if the genome is available, a simulated run can be created to resemble the results of amplicon sequencing with the specific panel required by the pipeline. We provide a tool that facilitates the creation of such data.  
Simulated data can also be useful 


## Rationale for selection of resources
Pipelines are usually tailored to a particular sequences technology, such as whole-genome or amplicon sequencing. While Pf8 is based on (selective) Whole Genome Sequencing (sWGS) data, GenRE Mekong utilised the SPotMalaria panel for amplicon sequencing (download available at [https://www.malariagen.net/wp-content/uploads/2023/11/20200705-GenRe-04b-SpotMalaria-SupplementaryFile1.xlsx](https://www.malariagen.net/wp-content/uploads/2023/11/20200705-GenRe-04b-SpotMalaria-SupplementaryFile1.xlsx).
Some GenRE Mekong samples have been re-analysed with sWGS for Pf8, thus providing data for the same sample from two different technologies and two independent resistance phenotype calls.
We will use the intersect between Pf8 and GenRE Mekong as the gold-standard dataset.

In Pf8, drug-resistance phenotypes are inferred based on criteria described in [this PDF](https://pf8-release.cog.sanger.ac.uk/Pf8_resistance_classification.pdf)


#### Rationale for selecting the dataset
#### Recipe for creating the dataset

### Simulated data