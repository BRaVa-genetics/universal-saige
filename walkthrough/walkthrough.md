# Walkthrough

## Introduction

Universal-saige has been created to standardise and ease the usage of [SAIGE](https://github.com/saigegit/SAIGE) for [BRaVa](https://brava-genetics.github.io/BRaVa/) across a variety of different computing environments.
In this walkthrough we will demonstrate how to generate gene and variant associations for the BRaVa phenotype on chromosome 22 with a final section on sanity-checking results. 

## Support

If at any point you run into issues or have any questions please create an issue in this (public) repository. You can also email `barney.hill (at) ndph.ox.ac.uk`

## Requirements

### Data

- Genotype data, plink (optional), ideally used in place of exome data for step 0
- Exome data, vcf or plink. 
- Sample IDs, a list of sampleIDs to analyse.
- Annotation file
- BRaVa phenotype file (provided)

### Environment

The only env requirement for this walkthrough is access to a linux machine with either Docker or Singularity available. With Docker or Singularity we can run Wei Zhou's [SAIGE Docker container](https://hub.docker.com/r/wzhou88/saige)
giving guarantees that analyses across cohorts are equivalent and easily reproducible. 

## Setup

To run universal-saige we need to download plink and the SAIGE container. These steps are seperated out into `download_resources.sh`:

### Setup (if using Docker)
`bash download_resources.sh --saige-image --plink`
### Setup (if using Singularity)
`bash download_resources.sh --saige-image --plink --singularity`

## Step 0 

To start we must generate the sparse GRM (genetic relatedness matrix) and processed plink files for usage in variance ratio estimation during step 1. While this step may take several hours to run it only has to be executed once per biobank.
Step 0 supports (genotype data, plink format), (exome data, vcf format) and (exome data, plink format) as inputs although we reccomend the usage of (genotype, plink format) in order to reduce runtime and maximise the number of independent sites.
For this step we reccomend using a larger machine - most functions in this step are parallelised across CPU cores and will benefit from high RAM. 

To begin clone the latest version of universal-saige
```
git clone git@github.com:BRaVa-genetics/universal-saige.git
cd universal-saige
mkdir out in
```

For this walkthrough we will be running step 0 with genotype plinks.
NOTE: Docker and Singularity require all input files to be within one directory that must not contain any linked files (so no `ln -s` your input files into your dir).
Currently my directory looks like:

```
.
├── ...
├── 00_step0_VR_and_GRM.sh
├── out/
├── in/
│   ├── ukb_genotypes_chr*.bed   # genotype bed files
│   ├── ukb_genotypes_chr*.bim.  # genotype bim files
│   ├── ukb_genotypes_chr*.fam.  # genotype fam files
│   ├── sample_ids.txt           # --sampleIDs
```

And I run step 0 with the arguments:
```
bash 00_step0_VR_and_GRM.sh \
    --geneticDataDirectory in/ \
    --geneticDataFormat "plink" \
    --geneticDataType "genotype" \
    --outputPrefix out/walkthrough \
    --sampleIDs in/sample_ids.txt \
    --generate_plink_for_vr \
    --generate_GRM
```

This took X hours with X cores and M GB memory. Checking my out directory I can now see:

```
.
├── ...
├── out/
│   ├── walkthrough.plink_for_var_ratio.bed
│   ├── walkthrough.plink_for_var_ratio.bed
│   ├── walkthrough.plink_for_var_ratio.bed
│   ├── walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx
│   ├── walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt

```

## Step 1

In step 1 we will be estimating the variance ratios for the association tests in step 2 (to be performed once per phenotype). For this walkthrough we'll use the continuous trait height as an example.

`head phenoFile.txt`

| IID | Type_two_diabetes | age | assessment_centre | PC1 |
|---|---|---|---|---|
| 3421 | 1 | 68 | 42 | 0.013412 |
| 4567 | 0 | 51 | 66 | -0.200134 |

```
bash 01_step1_fitNULLGLMM.sh \
    -t binary \
    --genotypePlink out/walkthrough.plink_for_var_ratio \
    --phenoFile in/phenoFile.txt \
    --phenoCol "female_infertility_binary" \
    --covarColList "age,assessment_centre" \
    --categCovarColList "assessment_centre" \
    --sampleIDs in/sample_ids.txt \
    --sampleIDCol "IID" \
    --outputPrefix out/walkthrough \
    --isSingularity false \
    --sparseGRM out/walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMID out/walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
```

## Step 2

```
bash 02_step2_SPAtests_variant_and_gene.sh \
    --chr chr1 \
    --testType "variant" \
    --plink "in/exome" \
    --modelFile out/walkthrough \
    --varianceRatio out/walkthrough \
    --groupFile in/ \
    --outputPrefix $output_prefix \
    --isSingularity false \
    --subSampleFile in/subsample_file/* \
    --sparseGRM in/GRM/* \
    --sparseGRMID in/GRM_samples/*
```
