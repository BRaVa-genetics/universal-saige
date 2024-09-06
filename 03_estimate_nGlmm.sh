binary_phenos=""
cont_phenos=""
PHENO_FILE=""
COVAR_LIST=""
SPARSE_GRM_FILE=""
SPARSE_GRM_ID_FILE=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --binaryPhenos)
      binary_phenos="$2"
      shift # past argument
      shift # past value
      ;;
    --contPhenos)
      cont_phenos="$2"
      shift # past argument
      shift # past value
      ;;
    --phenoFile)
      PHENO_FILE="$2" # past argument
      shift # past value
      ;;
    --covarList)
      COVAR_LIST="$2"
      shift # past argument
      shift # past value
      ;;
    --sparseGRM)
      SPARSE_GRM_FILE="$2"
      shift # past argument
      shift # past value
      ;;
    --sparseGRMID)
      SPARSE_GRM_ID_FILE="$2"
      shift # past argument
      shift # past value
      ;;
    -h|--help)
      echo "usage: 03_estimate_nGlmm.sh
  required:
    --binaryPhenos: space separated list of binary phenotypes.
    --contPhenos: space separated list of continuous phenotypes.
    --phenoFile: filename of the phenotype file. This must be relative to, and contained within, the current working directory.
    --sparseGRM: filename of the sparseGRM .mtx file. This must be relative to, and contained within, the current working directory.
    --sparseGRMID: filename of the sparseGRM ID file. This must be relative to, and contained within, the current working directory.
      "
      shift # past argument
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters

WD=$( pwd )

docker run -i \
  -e HOME=${WD} \
  -e BINARY_PHENOS="$binary_phenos" \
  -e CONT_PHENOS="$cont_phenos" \
  -e PHENO_FILE="$PHENO_FILE" \
  -e COVAR_LIST="$COVAR_LIST" \
  -e SPARSE_GRM_FILE="$SPARSE_GRM_FILE" \
  -e SPARSE_GRM_ID_FILE="$SPARSE_GRM_ID_FILE" \
  -v ${WD}/:$HOME/ \
  -v /mnt/project/:/mnt/project/ \
  wzhou88/saige:1.3.4 /bin/bash << EOF

set -x  # Enable debugging output

sed -i '/setgeno/d' R/SAIGE_extractNeff.R

# Install the package
R CMD INSTALL .

# Define phenotype variables (binary and continuous)
echo "pheno,nglmm" > \$HOME/neff.csv

echo "Estimating neff for cont phenotypes: \$CONT_PHENOS"
for pheno in \$CONT_PHENOS; do
    echo "Estimating neff for \$pheno"
    echo "PHENO_FILE: \$PHENO_FILE"
    echo "COVAR_LIST: \$COVAR_LIST"
    echo "SPARSE_GRM_FILE: \$SPARSE_GRM_FILE"
    echo "SPARSE_GRM_ID_FILE: \$SPARSE_GRM_ID_FILE"
    
    # Check if Rscript is available
    which Rscript || echo "Rscript not found in PATH"
    
    # Check if the R script file exists
    ls -l extdata/extractNglmm.R || echo "extractNglmm.R not found"
    
    # Run the Rscript command with error checking
    Rscript extdata/extractNglmm.R \
        --phenoFile \$PHENO_FILE \
        --phenoCol \$pheno \
        --covarColList \$COVAR_LIST \
        --traitType 'quantitative' \
        --sparseGRMFile \$SPARSE_GRM_FILE \
        --sparseGRMSampleIDFile \$SPARSE_GRM_ID_FILE \
        --useSparseGRMtoFitNULL TRUE 2>&1 | tee rscript_output.log
    
    # Check if the Rscript command was successful
    if [ \$? -ne 0 ]; then
        echo "Rscript command failed. Check rscript_output.log for details."
        cat rscript_output.log
    else
        # Process the output only if Rscript was successful
        grep 'Nglmm' rscript_output.log | awk -v pheno_var="\$pheno" '{print pheno_var "," \$2}' >> \$HOME/neff.csv
    fi
done

echo "Estimating neff for binary phenotypes: \$BINARY_PHENOS"
for pheno in \$BINARY_PHENOS; do
    echo "Estimating neff for \$pheno"
    echo "PHENO_FILE: \$PHENO_FILE"
    echo "COVAR_LIST: \$COVAR_LIST"
    echo "SPARSE_GRM_FILE: \$SPARSE_GRM_FILE"
    echo "SPARSE_GRM_ID_FILE: \$SPARSE_GRM_ID_FILE"
    
    # Check if Rscript is available
    which Rscript || echo "Rscript not found in PATH"

    # Check if the R script file exists
    ls -l extdata/extractNglmm.R || echo "extractNglmm.R not found"

    # Run the Rscript command with error checking
    Rscript extdata/extractNglmm.R \
        --phenoFile \$PHENO_FILE \
        --phenoCol \$pheno \
        --covarColList \$COVAR_LIST \
        --traitType 'binary' \
        --sparseGRMFile \$SPARSE_GRM_FILE \
        --sparseGRMSampleIDFile \$SPARSE_GRM_ID_FILE \
        --useSparseGRMtoFitNULL TRUE 2>&1 | tee rscript_output.log
    
    # Check if the Rscript command was successful
    if [ \$? -ne 0 ]; then
        echo "Rscript command failed. Check rscript_output.log for details."
        cat rscript_output.log
    else  
        # Process the output only if Rscript was successful
        grep 'Nglmm' rscript_output.log | awk -v pheno_var="\$pheno" '{print pheno_var "," \$2}' >> \$HOME/neff.csv
    fi
done

echo "Finished estimating neff"
echo "Contents of neff.csv:"
cat \$HOME/neff.csv

EOF