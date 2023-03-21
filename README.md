# Universal SAIGE

This code runs sets up SAIGE wzhou88/saige:1.1.6.3 docker or singularity container on your computer (if it is not already present). 

Note that your cluster or computer will require an internet connection for this. Also, either singularity or docker will need to be installed on your computer in order to grab the container.

When executing `main.sh`, you need to pass the following collection of filepaths and parameters.

phenoCol
subSampleFile
covarColList
qCovarColList
chrom # Phenotype+covariate variables

GMMATmodelFile
SAIGEOutputFile
sparseGRMFile
outputPrefix # Where to write files

exomeFileBim
exomeFileBed
exomeFileFam
vcfFile # Exome data
 
genotypeFileBim
genotypeFileBed
genotypeFileFam # Genotype data

The code then runs a series of sanity checks (see below in `checks.sh`).

## Cohort variation
- VCF or plink exome files
- Genotype or no genotype file for GRM
- Singularity or docker availability
- 
