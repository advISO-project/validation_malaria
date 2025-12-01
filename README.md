# A resource for malaria gold-standard validation data

## Aim of this resource
The aim of this resource is to provide datasets that may be used in validating pipelines for malaria genomics pipelines in order to provide evidence for accreditation of the pipeline and process.  
Two datasets are provided:  
- high-quality real-world data from public resources
- simulated data: example data and recipes

## Assumptions
Analysis pipelines differ in scope and methodology. Thus, no single dataset will be applicable to every possible genomics pipeline. For the selection of real-world public datasets, we have used the following assumptions:  

## Details about the data sets
This section provides details about the provided datasets, including recipes for creating the data, which may serve as a useful starting point for creating tailored datasets for different purposes.  

### Real-world data
The data in this section was obtained from two high-profile public malaria genomics resources:  
1. [MalariaGEN Pf8](https://www.malariagen.net/data_package/open-dataset-plasmodium-falciparum-v80/)
2. [GenRE Mekong](https://www.malariagen.net/resource/29/)

## Rationale for selection of resources
Pipelines are usually tailored to a particular sequences technology, such as whole-genome or amplicon sequencing. While Pf8 is based on (selective) Whole Genome Sequencing (sWGS) data, GenRE Mekong utilised the SPotMalaria panel for amplicon sequencing (download available at [https://www.malariagen.net/wp-content/uploads/2023/11/20200705-GenRe-04b-SpotMalaria-SupplementaryFile1.xlsx](https://www.malariagen.net/wp-content/uploads/2023/11/20200705-GenRe-04b-SpotMalaria-SupplementaryFile1.xlsx).
Some GenRE Mekong samples have been re-analysed with sWGS for Pf8, thus providing data for the same sample from two different technologies and two independent resistance phenotype calls.
We will use the intersect between Pf8 and GenRE Mekong as the gold-standard dataset.

In Pf8, drug-resistance phenotypes are inferred based on criteria described in [this PDF](https://pf8-release.cog.sanger.ac.uk/Pf8_resistance_classification.pdf)


#### Rationale for selecting the dataset
#### Recipe for creating the dataset

### Simulated data