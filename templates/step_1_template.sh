# Note that because of how the containers are set up, the files passed must be
# contained within the working directory. To change this, edit the
# `run_container.sh` script which calls the docker/singularity container for
# all steps.

TRAIT_TYPE={"quantitative","binary"}

# Path to the plink files created by step 0 (used to generated the variance 
# ratio). It will be of the form ${out_step0}.plink_for_var_ratio.{bim,bed,fam}
# where ${out_step0} is the output prefix to data generated by step 0. 
plink_for_var_ratio=["path/to/varianceratio/plink/files"]

# Path to the phenotype file. This file should be tab delimited, contain
# column names, and a column of sample IDs that matches the sample IDs of the
# genetic data. 
# Note that the phenotype file must also contain all of the covariates.
# Make sure not to include spaces or strange characters in the column names.
pheno_file=["path/to/phenotype/file"]

# The column name for the phenotype to be analysed - this must match exactly!
pheno=["pheno_col"]

# Comma separated (without spaces!) list of covariates. These covariate names
# must exactly match the corresponding column names in ${pheno_file}.
covariates=["age,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"]
# Comma separated (without spaces!) list of categorical covariates
# If a categorical variable is to be included as a covariate, it should be
# included in both ${covariates} and ${categorical_covariates}. 
categorical_covariates=["sex"]

# Optional! Path to a sample ID file. This should be a single column of sample
# IDs, with no header. The code will simply restrict the phenotype file to 
# these samples.
sample_id_path=["path/to/sampleIDs"]

# The column name for containing the sample IDs - this must match exactly!
# Note: the default is IID.
sample_id_col=["IID"]

# Output prefix for the data generated by this step (step 1).
out_step1=["out"]

# Is singularity available? If not, it is assumed that docker is available.
IS_SINGULARITY={"false","true"}

# Path to the GRM file created by step 0.
# It will be of the form 
# ${out_step0}.relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx where
# ${out_step0} is the output prefix to data generated by step 0. 
GRM=["path/to/GRM/file"]

# Path to the GRM sample IDs file created by step 0.
# It will be of the form 
# ${out_step0}.relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
# where ${out_step0} is the output prefix to data generated by step 0.
GRM_samples=["path/to/GRM/sampleID/file"]

bash 01_step1_fitNULLGLMM.sh \
    --t ${TRAIT_TYPE} \
    --genotypePlink ${plink_for_var_ratio} \
    --phenoFile ${pheno_file} \
    --phenoCol ${pheno} \
    --covarColList ${covariates} \
    --categCovarColList ${categorical_covariates} \
    --sampleIDs ${sample_id_path} \
    --sampleIDCol ${sample_id_col} \
    --outputPrefix ${out_step1} \
    --isSingularity ${IS_SINGULARITY} \
    --sparseGRM ${GRM} \
    --sparseGRMID ${GRM_samples}

