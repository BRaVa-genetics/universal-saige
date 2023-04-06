# Universal SAIGE

Run SAIGE preprocessing, step 1 and 2 without any hassle.

## System Requirements
- Internet connection (to download [SAIGE docker image](https://hub.docker.com/r/wzhou88/saige))
- Docker OR Singularity
- Linux OR Mac

## Input data (required)
- VCF or Plink WES data
- File containing sample ID's with the corresponding phenotype and covariate values.

## Input data (Optional)
- Genotyped sequences for every sample included in the WES data.

## Usage
### Step 0 (once per cohort/biobank)
Taking {vcf/plink/bgen} {genotype/WES/WGS} files generate plink file for the variance ratios and the GRM.  

```
usage: 00_step0_VR_and_GRM.sh
  required:
    --geneticDataDirectory: directory containing the genetic data (genotype/WES/WGS data in the format plink/vcf/bgen)
    --geneticDataFormat: format of the genetic data {plink,vcf,bgen}.
  optional:
    -o,--outputPrefix:  output prefix of the SAIGE step 1 output.
    -s,--isSingularity (default: false): is singularity available? If not, it is assumed that docker is available.
    --sampleIDCol (default: IID): column containing the sample IDs in the phenotype file, which must match the sample IDs in the plink files.
    --generate_GRM (default: false): generate GRM for the genetic data.
    --generate_plink_for_vr (default: false): generate plink file for vr.
```

### Step 1 (once per phenotype)

```
usage: 01_step1_fitNULLGLMM.sh
  required:
    -t,--traitType: type of the trait {quantitative,binary}.
    --genotypePlink: plink filename prefix of bim/bed/fam files. These must be present in the working directory at ./in/plink_for_vr_bed/
    --sparseGRM: filename of the sparseGRM .mtx file. This must be present in the working directory at ./in/sparse_grm/
    --sparseGRMID: filename of the sparseGRM ID file. This must be present in the working directory at ./in/sparse_grm/
    --phenoFile: filename of the phenotype file. This must be present in the working directory at ./in/pheno_files/
    --phenoCol: the column names of the phenotype to be analysed in the file specified in --phenoFile.
  optional:
    -o,--outputPrefix:  output prefix of the SAIGE step 1 output.
    -s,--isSingularity (default: false): is singularity available? If not, it is assumed that docker is available.
    -c,--covarColList: comma separated column names (e.g. age,pc1,pc2) of continuous covariates to include as fixed effects in the file specified in --phenoFile.
    --categCovarColList: comma separated column names of categorical variables to include as fixed effects in the file specified in --phenoFile.
    --sampleIDCol (default: IID): column containing the sample IDs in the phenotype file, which must match the sample IDs in the plink files.
```

### Step 2 (once per chromosome per phenotype)

```

usage: 02_step2_SPAtests_variant_and_gene.sh
  required:
    --testType: type of test {variant,group}.
    -p,--plink: plink filename prefix of bim/bed/fam files. These must be present in the working directory at ./in/plink_for_vr_bed/
    --vcf vcf exome file. If a plink exome file is not available then this vcf file will be used. These must be present in the working directory at ./in/vcf/
    --modelFile: filename of the model file output from step 1. This must be in relation to the working directory.
    --varianceRatio: filename of the varianceRatio file output from step 1. This must be in relation to the working directory.
    --sparseGRM: filename of the sparseGRM .mtx file. This must be present in the working directory at ./in/sparse_grm/
    --sparseGRMID: filename of the sparseGRM ID file. This must be present in the working directory at ./in/sparse_grm/
	--chr: chromosome to test.
  optional:
    -o,--outputPrefix:  output prefix of the SAIGE step 2 output.
    -s,--isSingularity (default: false): is singularity available? If not, it is assumed that docker is available.
    -g,--groupFile: required if group test is selected. Filename of the annotation file used for group tests. This must be in relation to the working directory.
```
