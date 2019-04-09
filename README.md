# Microbiome Mini Project
Project investigating the use of the R package dada2 for analysis of Illumina 16S RNA sequencing data.  Work includes going through the tutorial and extending use to dataset of dogs with inflammatory bowel disease vs healthy controls.  

## Table of Contents
1) Project Motivation <br>
2) File Description <br>
3) Libraries Required <br>
4) Summary of Results <br>
5) Licensing and Acknowledgements <br>

## Project Motivation
To learn how to analyze microbiome data in R.

## File Descriptions
Data files used for the dada2 tutorial can be found here:
https://benjjneb.github.io/dada2/tutorial.html

Data files for the canine IBD analysis can be found here:
https://www.ebi.ac.uk/ena/data/view/PRJEB13362

Files required for taxonomical assignment can be found here:
https://benjjneb.github.io/dada2/training.html

**README.md** - This file, describing the contents of this repo

**.gitignore** - The gitignore file

**dada2_tutorial.R** - R script detailing steps of the dada2 tutorial

**dada2_canine_ibd.R** - R script detailing the analysis of canine ibd files using the dada2 package

**results_summary.pdf** - Summary file showing graphs of alpha diversity and differential proportions of top 20 taxa

## Libraries Required
**dada2** - For accurate, high-resolution sample inference from amplicon sequencing data

**DECIPHER** - Alternative method for taxonomic classification

**phyloseq** and **ggplot2** - Handling, analysis and visualization of high-throughput microbiome census data

## Summary of results
See results_summary.pdf.  Dogs with IBD have reduced alpha diversity compared to healthy dogs.  Various differences in microbiome were detected between IBD and healthy dogs.

## Licenses and Acknowledgements
The paper describing the original analysis of the data can be found here: https://www.nature.com/articles/nmicrobiol2016177.  The GitHub profile of the original author's work can be found here: https://github.com/ElDeveloper/dogs/.  
