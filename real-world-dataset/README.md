# Pf8/GenRe Mekong gold-standard datasets for pipeline validation
This folder contains data files, created by this Jupyter notebook:  
[scripts/pf8genre.ipynb](scripts/pf8genre.ipynb)

## Aim
The aim of the notebook is to create a gold-standard dataset for malaria genomics pipelines that predict drug resistance.

The following resources will be used:
1. [MalariaGEN Pf8](https://www.malariagen.net/data_package/open-dataset-plasmodium-falciparum-v80/)
2. [GenRE Mekong](https://www.malariagen.net/resource/29/)

## Rationale for selection of resources
Pipelines are usually tailored to a particular sequencing  technology, such as whole-genome or amplicon sequencing. While Pf8 is based on (selective) Whole Genome Sequencing (sWGS) data, GenRE Mekong utilised the SPotMalaria panel for amplicon sequencing (download available at [https://www.malariagen.net/wp-content/uploads/2023/11/20200705-GenRe-04b-SpotMalaria-SupplementaryFile1.xlsx](https://www.malariagen.net/wp-content/uploads/2023/11/20200705-GenRe-04b-SpotMalaria-SupplementaryFile1.xlsx).
Some GenRE Mekong samples have been re-analysed with sWGS for Pf8, thus providing data for the same sample from two different technologies and two independent resistance phenotype calls.
We will use the intersect between Pf8 and GenRE Mekong as the gold-standard dataset.

In Pf8, drug-resistance phenotypes are inferred based on criteria described in [this PDF](https://pf8-release.cog.sanger.ac.uk/Pf8_resistance_classification.pdf)

## Executive Summary and results files
Three datasets of samples are created in this notebook. An overview of the files can be found in the [Summary-of-files-created section of the notebook](scripts/pf8genre.ipynb#Summary-of-files-created)