# Walkthrough

## Introduction

Universal-saige has been created to standardise and ease the usage of [SAIGE](https://github.com/saigegit/SAIGE) for [BRaVa](https://brava-genetics.github.io/BRaVa/) across a variety of different computing environments.
In this walkthrough we will demonstrate how to generate gene and variant associations for the BRaVa phenotype on chromosome 22 with a final section on sanity-checking results. 

## Support

If at any point you run into issues or have any questions please create an issue in this (public) repository. You can also email `barney.hill (at) ndph.ox.ac.uk`

## Requirements

### Data

- 

### Environment

The only env requirement for this walkthrough is access to a linux machine with either Docker or Singularity available. With Docker or Singularity we can run Wei Zhou's [SAIGE Docker container](https://hub.docker.com/r/wzhou88/saige)
giving guarantees that analyses across cohorts are equivalent and easily reproducible. 

## Step 0 

To begin we must generate the sparse GRM (genetic relatedness matrix) and processed plink files for usage in variance ratio estimation during step 1. While this step may take several hours to run it only has to be executed once per biobank.
Step 0 supports (genotype data, plink format), (exome data, vcf format) and (exome data, plink format) as inputs although we reccomend the usage of (genotype, plink format) in order to reduce runtime and maximise the number of independent sites.
For this step we reccomend using a larger machine - most functions in this step are parallelised across CPU cores and will benefit from high RAM. 

To begin we clone V1 of universal-saige
```
git clone git@github.com:BRaVa-genetics/universal-saige.git
cd universal-saige
```

```
bash 00_step0_VR_and_GRM.sh \
    --geneticDataDirectory "in/" \
    --geneticDataFormat "plink" \
    --geneticDataType "genotype" \
    --outputPrefix $out \
    --sampleIDs ~/in/sample_ids/* \
    --generate_plink_for_vr \
    --generate_GRM
```

## Step 1

```
bash 01_step1_fitNULLGLMM.sh \
    -t binary \
    --genotypePlink "genotype" \
    --phenoFile in/pheno_list/* \
    --phenoCol "female_infertility_binary" \
    --covarColList "age,assessment_centre" \
    --categCovarColList "assessment_centre" \
    --sampleIDs in/sample_ids/* \
    --sampleIDCol "IID" \
    --outputPrefix ${output_prefix} \
    --isSingularity false \
    --sparseGRM GRM \
    --sparseGRMID in/GRM_samples/*
```

## Step 2

```
bash 02_step2_SPAtests_variant_and_gene.sh \
    --chr $chrom \
    --testType "variant" \
    --plink "in/exome" \
    --modelFile in/model_file/* \
    --varianceRatio in/variance_ratio/* \
    --groupFile in/group_file/* \
    --outputPrefix $output_prefix \
    --isSingularity false \
    --subSampleFile in/subsample_file/* \
    --sparseGRM in/GRM/* \
    --sparseGRMID in/GRM_samples/*
```
