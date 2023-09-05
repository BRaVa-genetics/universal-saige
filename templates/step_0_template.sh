# Note that because of how the containers are set up, the files passed must be
# contained within the working directory. To change this, edit the
# `run_container.sh` script which calls the docker/singularity container for
# all steps.

# All files contained within --geneticDataDirectory of the type flagged by
# --geneticDataFormat will be globbed, so please ensure that this contains all
# of the autosomes for just one biobank/cohort and not multiple!
genetic_data_directory=["in/"]

# Note, VCF files must be gzipped with `.vcf.gz` file extensions.
GENETIC_DATA_FORMAT={"plink","vcf"}

GENETIC_DATA_TYPE={"WES","WGS","genotype"}

# Is singularity available? If not, it is assumed that docker is available.
IS_SINGULARITY={"false","true"}

# Output prefix for the data generated by this step (step 0).
out_step0=["out"]

# These are the sampleIDs (single column) to be used to define the GRM.
# Note that if nothing is passed (the string is empty), then all of the samples
# in the files will be used
sample_id_path=["path/to/sampleIDs"]

# The inclusion of the following two flags ensures that the variance ratio file
# and GRM files are created.
# --generate_plink_for_vr (used in step 1 and flagged by --genotypePlink).
# --generate_GRM (used in step 1 and flagged by --sparseGRM and --sparseGRMID).

bash 00_step0_VR_and_GRM.sh \
    --geneticDataDirectory ${genetic_data_directory} \
    --geneticDataFormat ${GENETIC_DATA_FORMAT} \
    --geneticDataType ${GENETIC_DATA_TYPE} \
    --isSingularity ${IS_SINGULARITY} \
    --outputPrefix ${out_step0} \
    --sampleIDs ${sample_id_path} \
    --generate_plink_for_vr \
    --generate_GRM
