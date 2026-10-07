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
  * [Setup (if using Singularity or Apptainer)](#setup-if-using-singularity-or-apptainer)
* [Step 0](#step-0)
* [Step 1](#step-1)
* [Step 2](#step-2)
* [Step 3](#step-3)

## Introduction

Universal-SAIGE has been created to standardise the usage of [SAIGE](https://github.com/saigegit/SAIGE) for [BRaVa](https://brava-genetics.github.io/BRaVa/) across a variety of different computing environments.

In this walkthrough we will demonstrate how to generate gene and variant associations for the BRaVa phenotype HDL cholesterol on chromosome 11 with a final section on sanity-checking results. 

## Support

If at any point you run into issues or have any questions please create an issue in this (public) repository. You can also email `barney.hill@ndph.ox.ac.uk` or `duncan.palmer@ndph.ox.ac.uk`.

## Caution
> [!WARNING]
> A few things to be aware of:
> - A binary trait with fewer than 100 cases (in the fitted samples), a categorical covariate level with fewer than 10 cases or controls, a separated covariate model, or a fit that did not converge is refused in step 1, and again in step 2, with the reason in the log. That refusal is the right answer for such a trait; `SAIGE_FIT_GATES=0` turns it into a warning (not recommended).
> - When running sex-specific phenotypes do not include sex as a covariate. This can cause invalid results/crashes. Drop every sex term (All of Us used `age,age2,PC1,...,PC20` and no categorical covariate), and pass `--sex F` or `--sex M` to step 1, which then fits only that sex and drops the rest (SAIGE's `--FemaleOnly`/`--MaleOnly`, as All of Us did). A numeric `sex` column is read as 0 = female, 1 = male, BRaVa's coding; add `--femaleCode 1 --maleCode 0` if yours is the other way round, and check the counts step 1 prints. Where the genotype file has genetic sex (`.fam` column 5, or the `.psam` `SEX` column; PLINK fixes 1 = male, 2 = female), step 1 checks the phenotype file's sex against it: more than half disagreeing means the codes are flipped and is refused, more than 1% is warned loudly, and a genotype file without usable sex (e.g. converted from VCF) is warned loudly as unverifiable. The phenotype file's `sex` always decides who is fitted.

## Requirements

### Data

- Genotype data, plink (optional), ideally used in place of exome data for step 0
- Exome data in PLINK 2 (`.pgen/.pvar/.psam`, recommended) or PLINK 1 (`.bed/.bim/.fam`) format; a VCF is converted once with plink2 (`plink2 --vcf exome.chr11.vcf.gz --make-pgen --out exome.chr11`). 
- Sample IDs, (ancestry specific)
- Annotation (group) file, generated [here](https://github.com/BRaVa-genetics/variant-annotation) ([details](https://docs.google.com/document/d/1emWqbX8ohi-9rYIW_pKSAFiMHZZUV6zyXwg7qWJNdlc/edit#heading=h.puz6ua3vxnca); [thresholds and versions](https://docs.google.com/document/d/11Nnb_nUjHnqKCkIB3SQAbR6fl66ICdeA-x_HyGWsBXM/edit#heading=h.649be2dis6c1))
- BRaVa phenotype file (.tsv) with 'IID' (sample ID) column and covariates

### Environment

The only env requirement for this walkthrough is access to a linux machine with Docker, Singularity or Apptainer available (Apptainer is Singularity's successor; `--isSingularity true` uses whichever of `singularity` and `apptainer` is on PATH). With any of them we run the slim SAIGE build (`astheeggeggs/saige-slim` on Docker Hub, pinned by tag in `download_resources.sh`), which gives the same guarantee that analyses across cohorts are equivalent and reproducible. 

## Setup
To run universal-saige we need to download plink and the SAIGE image. These steps are separated out into `download_resources.sh`:

### Setup (if using Docker)
```
bash download_resources.sh --saige-image --plink2 --plink
```
### Setup (if using Singularity or Apptainer)
```
bash download_resources.sh --saige-image --plink2 --plink --singularity
```
The image is pulled without a cache and unpacked next to `resources/saige.sif`, not in `/tmp` or `$HOME`, which are often small on clusters. Set `APPTAINER_TMPDIR` to unpack it somewhere else.

## Step 0 
To start we must generate the sparse genetic relatedness matrix (GRM) and processed plink files for usage in variance ratio estimation during step 1. While this step may take several hours to run, it only has to be executed once per biobank/cohort.

Step 0 takes genotype or exome data in PLINK 1 (`--geneticDataFormat plink`) or PLINK 2 (`--geneticDataFormat pgen`) format, as steps 1 and 2 do; convert a VCF once with plink2 first. We recommend genotype array data here, to reduce runtime and maximise the number of independent sites.

For this step we recommend using a larger machine - most functions in this step are parallelised across CPU cores and will benefit from high RAM. 

To begin, clone the latest version of universal-saige
```
git clone git@github.com:BRaVa-genetics/universal-saige.git
cd universal-saige
mkdir -p out in/genotypes
```

For this walkthrough we will be running step 0 with plink files based on genotype array data, in a directory of their own: step 0 merges every `.bed` (or `.pgen`) in `--geneticDataDirectory`, so the exome files used in step 2 must not be in it. sample_ids.txt is a file with newline separated sample IDs.

> [!NOTE]
> Docker and Singularity require all input files to be within one directory that must not contain any linked files (so no `ln -s` your input files into your dir).

Currently my directory looks like:

```
.
├── ...
├── 00_step0_VR_and_GRM.sh
├── out/
├── in/
│   ├── genotypes/
│   │   ├── ukb_genotypes_chr*.bed   # genotype bed files
│   │   ├── ukb_genotypes_chr*.bim   # genotype bim files
│   │   ├── ukb_genotypes_chr*.fam   # genotype fam files
│   ├── sample_ids.txt               # --sampleIDs
```

And I run step 0 with the arguments:
```
bash 00_step0_VR_and_GRM.sh \
    --geneticDataDirectory in/genotypes/ \
    --geneticDataFormat "plink" \
    --geneticDataType "genotype" \
    --outputPrefix out/walkthrough \
    --sampleIDs in/sample_ids.txt \
    --generate_plink_for_vr \
    --generate_GRM
```

The GRM keeps pairs related at 0.05 or more (`--relatednessCutoff`, default 0.05), and its file name records the value. The same value must be passed to steps 0, 1 and 2; nothing in SAIGE checks that they agree. All of Us used 0.05, and 0.125 for its admixed amr cohort, whose GRM was too dense to fit at 0.05. Steps 0 and 1 print the GRM's mean number of relatives per sample (from the file header, so instantly), warn loudly above 100 (step 1 then refuses the GRM unless it is passed `--forceDenseGRM`), and warn when the GRM was built at a different cutoff from the step's; above 100 a fit can run for hours or never finish (All of Us amr: ~644 at 0.05, ~3.9 at 0.125; its other cohorts ~0.6).

This took 5 hours with 64 cores and 512 GB memory (for ~400K samples). Inspecting the `out/` directory, we can see:
```
.
├── ...
├── out/
│   ├── walkthrough.plink_for_var_ratio.bed
│   ├── walkthrough.plink_for_var_ratio.bim
│   ├── walkthrough.plink_for_var_ratio.fam
│   ├── walkthrough.plink_for_grm.{bed,bim,fam}   # the LD-pruned markers the GRM was built from
│   ├── walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx
│   ├── walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
```

The last line step 0 prints is the GRM's density, e.g. `sparse GRM: 400000 samples, ... relatives per sample (built at --relatednessCutoff 0.05)`; a loud warning there means the cutoff should be raised before step 1.

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
    --covarColList "age,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20" \
    --categCovarColList "sex" \
    --sampleIDs in/sample_ids.txt \
    --sampleIDCol "IID" \
    --outputPrefix out/HDL_cholesterol \
    --isSingularity false \
    --sparseGRM out/walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMID out/walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
```
> [!WARNING]
> A few things to note here:
> - The column names flagged in `--phenoCol`, `--covarColList` and `--categCovarColList` must _exactly_ match the column names in the filepath flagged by `--phenoFile`
> - The comma separated list of covariates flagged by `--covarColList` and `--categCovarColList` should not contain spaces (e.g. `age,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20`)
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

The labels in the `anno` lines are what `--annotations` refers to. Step 2 checks that every label you ask for appears in the group file, and refuses the run (listing the labels the file does have) if one does not. A group file written with other labels needs `--annotations` changed to match.

The top of the file looks like this:

`head in/ukb_brava_annotations.txt`

```
ENSG00000187634 var chr1:943315:T:C chr1:962890:T:A
ENSG00000187634 anno damaging_missense_or_protein_altering non_coding
ENSG00000187961 var chr1:961514:T:C chr1:962037:C:T chr1:962807:T:C 
ENSG00000187961 anno synonymous damaging_missense_or_protein_altering pLoF
```

Here, each gene (coded according to ensembl ID in column 1) receives two lines, a variant line (`var`) and an annotation line `anno` (column two). All subsequent information on each pair of gene specific lines contains space delimited information mapping the variant information onto the associated annotation(s). 

Finally, we perform the association testing for chromosome 11:

```
bash 02_step2_SPAtests_variant_and_gene.sh \
    --chr chr11 \
    --testType "group" \
    --plink in/ukb_wes_450k.qced.chr11 \
    --modelFile out/HDL_cholesterol.rda \
    --varianceRatio out/HDL_cholesterol.varianceRatio.txt \
    --groupFile in/ukb_brava_annotations.txt \
    --outputPrefix out/chr11_HDL_cholesterol \
    --annotations "pLoF,damaging_missense_or_protein_altering,other_missense_or_protein_altering,synonymous,pLoF:damaging_missense_or_protein_altering,pLoF:damaging_missense_or_protein_altering:other_missense_or_protein_altering:synonymous" \
    --sparseGRM out/walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMID out/walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
```
FlexRV on the same chromosome, once the AlphaMissense weights have been added to the group file
(`bash download_resources.sh --alphamissense`, then `bash 04_flexrv_groupfile.sh --group in/ukb_brava_annotations.txt --chr 11 --name AM --out in/ukb_brava_annotations.flexrv_AM.txt`):
```
bash 02_step2_SPAtests_variant_and_gene.sh \
    --chr chr11 \
    --plink in/ukb_wes_450k.qced.chr11 \
    --modelFile out/HDL_cholesterol.rda \
    --varianceRatio out/HDL_cholesterol.varianceRatio.txt \
    --groupFile in/ukb_brava_annotations.flexrv_AM.txt \
    --flexRVscore AM \
    --outputPrefix out/chr11_HDL_cholesterol.flexrv_AM \
    --sparseGRM out/walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMID out/walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
```
The pooled FlexRV p per gene is the `Group == Cauchy` row's `Pvalue_Burden`; the other rows are the transform sets.

> [!WARNING]
> The chromosome name flagged by `--chr` must _exactly_ match the chromosome name in the first column of the `.bim` (or `.pvar`): if it is labelled '11', pass `--chr 11`, not `--chr chr11`. Step 2 refuses a mismatch before SAIGE starts, and the message shows the file's own spelling.

This command took 1 hour 47 minutes (on an earlier SAIGE build). Step 2 runs single-threaded, so run chromosomes, and phenotypes, side by side rather than giving one run more cores. For verification of rare variant association results [genebass](https://app.genebass.org/) is a useful resource. Checking [HDL cholesterol](https://app.genebass.org/gene/undefined/phenotype/continuous-30760-both_sexes--irnt?resultIndex=gene-manhattan&resultLayout=full) we can see that APOC3 (ENSG00000110245) (pLoF, SKAT-O) has a association with $P=1.24\times 10^{-322}$. Looking at the gene result file `out/chr11_HDL_cholesterol.txt` we see the result:

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

    def get_lambda_gc(p):
        # median 1-df chi-square of the observed p-values over its expectation under the null
        return np.median(stats.chi2.isf(p, df=1))/stats.chi2.ppf(q=0.5, df=1)

    if type == "gene":
        pvals = results["Pvalue"][
                                  (results["max_MAF"] == max_maf) & 
                                  (results["Group"] == anno)]

    elif type == "variant":
        pvals = results["p.value"]

    if len(pvals) == 0:
        print(f"No results for {pheno}")
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

qqplot(gene_results, "HDL_cholesterol", "gene", max_maf=0.01, anno="damaging_missense_or_protein_altering")

variant_results = results_dir + "chr11_HDL_cholesterol.txt.singleAssoc.txt"
variant_results = pd.read_csv(variant_results, sep="\t")

qqplot(variant_results, "HDL_cholesterol", "variant")
```

The drivers run with fastTest off (the All of Us choice), so every _P_-value is computed in full, including those above 0.05, and $\lambda_{GC}$ can be read off the whole distribution.

<img src="https://user-images.githubusercontent.com/43707014/236252715-93df0a07-9799-4e50-85af-c679631a4bc3.png" width="500">

Taking a closer look:

`qqplot(variant_results[variant_results["p.value"] > 5E-8], "HDL_cholesterol", "variant")`

<img src="https://user-images.githubusercontent.com/43707014/236253166-f298e828-1954-4edf-96c9-c7638032dde9.png" width="500">

In this QQ-plot while we see some inflation from the expected p-values this is plausibly polygenicity given what we know about the trait. 

## Step 3

Finally, the effective sample size of each phenotype's null model, Nglmm, which the BRaVa meta-analysis uses. It is computed once per phenotype from the sparse GRM, on the samples step 1 fitted:

```
bash 03_estimate_nGlmm.sh \
    --contPhenos "HDL_cholesterol" \
    --phenoFile in/phenoFile.txt \
    --covarList "age,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20" \
    --sparseGRM out/walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMID out/walkthrough_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
    --outputFile out/neff.csv
```

`out/neff.csv` has one `pheno,nglmm` row per phenotype (`--binaryPhenos` takes the binary ones, space separated). Pass `--relatednessCutoff` if steps 0-2 used a value other than 0.05: Nglmm is computed at the same cutoff as the fit.
