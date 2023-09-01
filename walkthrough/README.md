# Universal-SAIGE walkthrough

### Contents
* [Introduction](#introduction)
* [Support](#support)
* [Caution](#caution)
* [Requirements](#requirements)
  * [Data](#data)
  * [Environment](#environment)
* [Setup](#setup)
  * [Setup (if using Docker)](#setup-if-using-docker)
  * [Setup (if using Singularity)](#setup-if-using-singularity)
* [Step 0](#step-0)
* [Step 1](#step-1)
* [Step 2](#step-2)

## Introduction

Universal-SAIGE has been created to standardise the usage of [SAIGE](https://github.com/saigegit/SAIGE) for [BRaVa](https://brava-genetics.github.io/BRaVa/) across a variety of different computing environments.

In this walkthrough we will demonstrate how to generate gene and variant associations for the BRaVa phenotype HDL cholesterol on chromosome 11 with a final section on sanity-checking results. 

## Support

If at any point you run into issues or have any questions please create an issue in this (public) repository. You can also email `barney.hill@ndph.ox.ac.uk` or `duncan.palmer@ndph.ox.ac.uk`.

## Caution
> **Warning**
> A few things to be aware of:
> - SAIGE can fail in a variety of ways due to low case count - we don't handle this within universal-SAIGE but step1/step2 failing across an entire phenotype x ancestry is a likely indicator for this
> - When running sex-specific phenotypes do not include sex as a covariate. This can cause invalid results/crashes

## Requirements

### Data

- Genotype data, plink (optional), ideally used in place of exome data for step 0
- Exome data, VCF or plink. 
- Sample IDs, (ancestry specific)
- Annotation file ([details found here](https://docs.google.com/document/d/1emWqbX8ohi-9rYIW_pKSAFiMHZZUV6zyXwg7qWJNdlc/edit#heading=h.puz6ua3vxnca](https://docs.google.com/document/d/11Nnb_nUjHnqKCkIB3SQAbR6fl66ICdeA-x_HyGWsBXM/edit#heading=h.649be2dis6c1)))
- BRaVa phenotype file (.tsv) with 'IID' (sample ID) column and covariates

### Environment

The only env requirement for this walkthrough is access to a linux machine with either Docker or Singularity available. With Docker or Singularity we can run Wei Zhou's [SAIGE Docker container](https://hub.docker.com/r/wzhou88/saige), giving guarantees that analyses across cohorts are equivalent and easily reproducible. 

## Setup
To run universal-saige we need to download plink and the SAIGE image. These steps are separated out into `download_resources.sh`:

### Setup (if using Docker)
```
bash download_resources.sh --saige-image --plink
```
### Setup (if using Singularity)
```
bash download_resources.sh --saige-image --plink --singularity
```

## Step 0 
To start we must generate the sparse genetic relatedness matrix (GRM) and processed plink files for usage in variance ratio estimation during step 1. While this step may take several hours to run, it only has to be executed once per biobank/cohort.

Step 0 supports (genotype data, plink format), (exome data, VCF format) and (exome data, plink format) as inputs although we reccomend the usage of (genotype, plink format) in order to reduce runtime and maximise the number of independent sites.

For this step we recommend using a larger machine - most functions in this step are parallelised across CPU cores and will benefit from high RAM. 

To begin, clone the latest version of universal-saige
```
git clone git@github.com:BRaVa-genetics/universal-saige.git
cd universal-saige
mkdir out in
```

For this walkthrough we will be running step 0 with plink files based on genotype array data. sample_ids.txt is a file with newline separated sample IDs.

> **Note**
> Docker and Singularity require all input files to be within one directory that must not contain any linked files (so no `ln -s` your input files into your dir).

Currently my directory looks like:

```
.
├── ...
├── 00_step0_VR_and_GRM.sh
├── out/
├── in/
│   ├── ukb_genotypes_chr*.bed   # genotype bed files
│   ├── ukb_genotypes_chr*.bim   # genotype bim files
│   ├── ukb_genotypes_chr*.fam   # genotype fam files
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

This took 5 hours with 64 cores and 512 GB memory (for ~400K samples). Inspecting the `out/` directory, we can see:
```
.
├── ...
├── out/
│   ├── walkthrough.plink_for_var_ratio.bed
│   ├── walkthrough.plink_for_var_ratio.bim
│   ├── walkthrough.plink_for_var_ratio.fam
│   ├── walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx
│   ├── walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt

```

## Step 1

In step 1 we will be fitting the null model for the association tests in step 2 (to be performed once per phenotype). For this walkthrough we'll use the continuous trait HDL cholesterol as an example. 

```
head in/phenoFile.txt
```

| IID | HDL_cholesterol | age | PC1       | ... |
| --- | --------------- | --- | --------- | --- |
| 3421 | 0.422          | 68  | 0.013412  | ... |
| 4567 | 0.342          | 51  | -0.200134 | ... |

```
bash 01_step1_fitNULLGLMM.sh \
    -t quantitative \
    --genotypePlink out/walkthrough.plink_for_var_ratio \
    --phenoFile in/phenoFile.txt \
    --phenoCol "HDL_cholesterol" \
    --covarColList "age,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" \
    --categCovarColList "sex" \
    --sampleIDs in/sample_ids.txt \
    --sampleIDCol "IID" \
    --outputPrefix out/HDL_cholesterol \
    --isSingularity false \
    --sparseGRM out/walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMID out/walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
```
> **Warning**
> A few things to note here:
> - The column names flagged in `--phenoCol`, `--covarColList` and `--categCovarColList` must _exactly_ match the column names in the filepath flagged by `--phenoFile`
> - The comma separated list of covariates flagged by `--covarColList` and `--categCovarColList` should not contain spaces (e.g. `age,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10`)
> - If a categorical variable is to be included as a covariate, it should be flagged by _both_ `--covarColList` and `--categCovarColList` (e.g. `sex` in the above command)
  
This command took 10 minutes with 4 cores. Checking the `out/` directory we can see:

```
.
├── ...
├── out/
│   ├── walkthrough.plink_for_var_ratio.bed
│   ├── walkthrough.plink_for_var_ratio.bim
│   ├── walkthrough.plink_for_var_ratio.fam
│   ├── walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx
│   ├── walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
│   ├── HDL_cholesterol.rda
│   ├── HDL_cholesterol.varianceRatio.txt
```

## Step 2

Step 2 requires variant annotations which can be generated [here](https://github.com/BRaVa-genetics/variant-annotation). A summary of the thresholds and software versioning used for variant annotation within BRaVa can be found [here](https://docs.google.com/document/d/11Nnb_nUjHnqKCkIB3SQAbR6fl66ICdeA-x_HyGWsBXM/edit#heading=h.649be2dis6c1), but you don't need to worry about the annoying version alignment if you follow our [steps](https://github.com/BRaVa-genetics/variant-annotation).

The top of the file looks like this:

`head in/ukb_brava_annotations.txt`

```
ENSG00000187634 var chr1:943315:T:C chr1:962890:T:A
ENSG00000187634 anno damaging_missense non_coding
ENSG00000187961 var chr1:961514:T:C chr1:962037:C:T chr1:962807:T:C 
ENSG00000187961 anno synonymous damaging_missense pLoF
```

Here, each gene (coded according to ensembl ID in column 1) receives two lines, a variant line (`var`) and an annotation line `anno` (column two). All subsequent information on each pair of gene specific lines contains space delimited information mapping the variant information onto the associated annotation(s). 

Finally, we perform the association testing for chromosome 11:

```
bash 02_step2_SPAtests_variant_and_gene.sh \
    --chr chr11 \
    --testType "group" \
    --plink in/ukb_wes_450k.qced.chr11.bed \
    --varianceRatio out/walkthrough \
    --groupFile out/walkthrough \
    --outputPrefix out/chr11_HDL_cholesterol \
    --annotations in/ukb_brava_annotations.txt \
    --isSingularity false \
    --subSampleFile in/sample_ids.txt \
    --sparseGRM out/walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMID out/walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
```
> **Warning**
> There's one more 'gotcha' here - you'll need to ensure that the chromosome name flagged by `--chr` _exactly_ matches the chromosome name in the .bim file. For example, if the chromosome is labelled as '11' in the first column of the .bim, `--chr chr11` will not work (but `--chr 11` will).

This command took 1 hour 47 minutes with 8 cores. For verification of rare variant association results [genebass](https://app.genebass.org/) is a useful resource. Checking [HDL cholesterol](https://app.genebass.org/gene/undefined/phenotype/continuous-30760-both_sexes--irnt?resultIndex=gene-manhattan&resultLayout=full) we can see that APOC3 (ENSG00000110245) (pLoF, SKAT-O) has a association with $P=1.24\times 10^{-322}$. Looking at the gene result file `out/chr11_HDL_cholesterol.txt` we see the result:

```
Region	Group	max_MAF	Pvalue	Pvalue_Burden	Pvalue_SKAT	BETA_Burden	SE_Burden	MAC	Number_rare	Number_ultra_rare
ENSG00000110245	pLoF	0.0100	3.318754e-305	4.741078e-306	5.078064e-288	0.034034	0.000910	1753.0	2.0	0.0
```
Replication! Of course we are using approximately the same cohort here (UK Biobank, European) but if you are following along with a HDL Cholesterol phenotype you will hopefully be able to observe similar results given sufficient power.

Another method of verification we reccomend is checking the QQ-plot, the expected vs observed _P_-values given the null hypothesis of the test. Below we plot the variant QQ-plot using Python:

```python
# Plot qqplot:
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import pandas as pd

def qqplot(results, pheno, type, max_maf=None, anno=None):

    def get_expected(n):
        exp = -np.log10(np.linspace(start=1,stop=1/n,num=n))
        return exp

    # Get 95% confidence interval
    def get_CI_intervals(n, CI=0.95):
        k = np.arange(1,n+1)
        a = k
        b = n+1-k
        intervals=stats.beta.interval(CI, a, b)
        return intervals

    def get_lambda_gc(chisq_vec):
        return np.median(chisq_vec)/stats.chi2.ppf(q=0.5, df=1)

    if type == "gene":
        pvals = results["Pvalue"][
                                  (results["max_MAF"] == max_maf) & 
                                  (results["Group"] == anno)]

    elif type == "variant":
        pvals = results["p.value"]

    if len(pvals) == 0:
        print(f"No results for {pheno} {sex}")
        return

    pvals = np.sort(pvals)
    pvals = pvals[pvals > 0]
    n = len(pvals)

    exp = get_expected(n)
    intervals = get_CI_intervals(n)

    x = exp[::-1]
    y = -np.log10(pvals)

    plt.figure(figsize=(10,10))
    if type == "gene": 
        plt.title(f"SAIGE-{type} results: {pheno} \nMax Allele Frequency:{max_maf} annotation:{anno} lambda_gc:{get_lambda_gc(pvals):.2f}")
    elif type == "variant": 
        plt.title(f"SAIGE-{type} results: {pheno}\nlambda_gc:{get_lambda_gc(pvals):.2f}")

    plt.xlabel("Expected -log10(p)")
    plt.ylabel("Observed -log10(p)")

    plt.fill_between(x=exp[::-1], y1=-np.log10(intervals[0]), y2=-np.log10(intervals[1]), color="gray", alpha=0.3, label="95% CI")

    plt.plot([0, max(exp[::-1])], [0, max(exp[::-1])], color="red")
    plt.plot(x, y, 'o')

results_dir = "out/"

gene_results = results_dir + "chr11_HDL_cholesterol.txt"
gene_results = pd.read_csv(gene_results, sep="\t")

qqplot(gene_results, "HDL_cholesterol", "gene", max_maf=0.01, anno="damaging_missense")

variant_results = results_dir + "chr11_HDL_cholesterol.txt.singleAssoc.txt"
variant_results = pd.read_csv(variant_results, sep="\t")

qqplot(variant_results, "HDL_cholesterol", "variant")
```

Note that due to fast testing enables results with $P > 0.05$ may be skewed and affect the $\lambda_{GC}$ value. 

<img src="https://user-images.githubusercontent.com/43707014/236252715-93df0a07-9799-4e50-85af-c679631a4bc3.png" width="500">

Taking a closer look:

`qqplot(variant_results[variant_results["p.value"] > 5E-8], "HDL_cholesterol", "variant")`

<img src="https://user-images.githubusercontent.com/43707014/236253166-f298e828-1954-4edf-96c9-c7638032dde9.png" width="500">

In this QQ-plot while we see some inflation from the expected p-values this is plausibly polygenicity given what we know about the trait. 
