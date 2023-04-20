TODO: 
- change step 1 to use the plink file for VR to restrict samples use awk to extract the right column and pass to SAIGE
- remove GRM and VR creation from step 1, and create an error message telling the user to run step 0.
- ensure that all samples in the vcf are present in the plink files
- code to generate the groupFiles
- step 2: add description of what the SUBSAMPLES file is and it's format to this document
- ensure step 0 help matches this help
- ensure step 1 help matches this help
- ensure step 2 help matches this help
- add a proposed clean directory structure for inputs and outputs

# Universal SAIGE

Run SAIGE preprocessing and steps 1 and 2 without any hassle.

## System Requirements
- Internet connection (to download [SAIGE docker image](https://hub.docker.com/r/wzhou88/saige))
- Docker OR Singularity
- Linux OR Mac

## Input data (required)
- `{WES, WGS}` data in plink (`.bim/.bed/.fam`) or VCF format
- Sample ID file containing the corresponding phenotype and covariate values.

## Input data (optional)
- Genotyping array data for every sample included in the WES data above.

## Usage
### Step 0 (once per cohort/biobank)
Take genotyping array data in `plink` format, or `{WES, WGS}` files in `{vcf, plink}` format, and generate variance ratios and a sparse GRM.

```
usage: 00_step0_VR_and_GRM.sh
```
required:
- `--geneticDataDirectory`: directory containing the genetic data (genotyping array data in `plink` format, or `{WES, WGS}` files in `{vcf, plink}` format)
- `--geneticDataFormat`: format of the genetic data `{vcf, plink}`.
- `--sampleIDs`: `.fam` file of the sample IDs that are present in the `{WES, WGS}` data. Note, if this is not _all_ of the samples in the `{WES, WGS}` dataset, the `{WES, WGS}` data must be filtered to these samples before running step 1.

optional:
- `-o`,`--outputPrefix`: output prefix from this program (SAIGE step 0) to be used as SAIGE step 1 input.
- `-s`,`--isSingularity` (default: `false`): is singularity available? If not, it is assumed that docker is available.
- `--generate_GRM` (default: false): generate GRM for the genetic data.
- `--generate_plink_for_vr` (default: false): generate plink file for vr.

### Step 1 (once per phenotype)

```
usage: 01_step1_fitNULLGLMM.sh
```
required:
- `-t`,`--traitType`: type of the trait `{quantitative, binary}`.
- `--genotypePlink`: variance ratio plink filename prefix of `.bim/.bed/.fam` files. This must relative to the current working directory. Note that samples will be restricted to samples present within the plink `.fam` file.
- `--sparseGRM`: filename of the sparseGRM `.mtx` file (output from step 0). This must be relative to the current working directory.
- `--sparseGRMID`: filename of the sparseGRM ID file (output from step 0). This must be relative to the current working directory.
- `--phenoFile`: filename of the phenotype file. This must be relative to the working directory.
- `--phenoCol`: the column names of the phenotype to be analysed in the file specified in `--phenoFile`.

optional:
- `-o`,`--outputPrefix`:  output prefix from this program (SAIGE step 1) to be used as SAIGE step 2 input.
- `-s`,`--isSingularity`: (default: false): is singularity available? If not, it is assumed that docker is available.
- `-c`,`--covarColList`: comma separated column names (e.g. `age,pc1,pc2`) of continuous covariates to include as fixed effects in the file specified in `--phenoFile`. Recall, proposed pilot fixed effect covariates are `age,age2,sex,age*sex,age2*sex,PCs`.
- `--categCovarColList`: comma separated column names of categorical variables to include as fixed effects in the file specified in --phenoFile.
- `--sampleIDCol` (default: IID): column containing the sample IDs in the phenotype file, which must match the sample IDs in the plink files.

### Step 2 (once per chromosome per phenotype)

```
usage: 02_step2_SPAtests_variant_and_gene.sh
```
required:
- `--testType`: type of test `{variant,group}`.
- `-p`,`--plink`: plink filename prefix of `.bim/.bed/.fam` for WES (or WGS restricted to exons). These must be relative to the current working directory.
- `--vcf` vcf exome file. If a set of plink files for the WES (or WGS restricted to exons) is not available then this vcf file will be used. This must be present in the current working directory.
- `--modelFile`: filename of the model file output from SAIGE step 1. This must be relative to the current working directory.
- `--varianceRatio`: filename of the varianceRatio file output from SAIGE step 1. This must be relative to the current working directory.
- `--sparseGRM`: filename of the sparseGRM `.mtx` file output from SAIGE step 0. This must be relative to the current working directory.
- `--sparseGRMID`: filename of the sparseGRM ID file output from SAIGE step 0. This must be relative to the current working directory.
- `--chr`: chromosome to test.

optional:
- `-o`,`--outputPrefix`: output prefix from this program (SAIGE step 2).
- `-s`,`--isSingularity` (default: false): is singularity available? If not, it is assumed that docker is available.
- `-g`,`--groupFile`: required if group test is selected. Filename of the annotation file used for group tests. This must be in relation to the working directory.
