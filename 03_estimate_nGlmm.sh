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


docker run -i -e HOME=${WD} -v ${WD}/:$HOME/ -v /mnt/project/:/mnt/project/ wzhou88/saige:1.3.4 /bin/bash << 'EOF'

sed -i '/setgeno/d' R/SAIGE_extractNeff.R

# Install the package
R CMD INSTALL .

# Define phenotype variables (binary and continuous)

echo "pheno,nglmm" >> neff.csv

echo "Estimating neff for binary phenotypes: $binary_phenos"
for pheno in $binary_phenos; do
    echo "Estimating neff for $pheno"
    Rscript extdata/extractNglmm.R \
        --phenoFile $PHENO_FILE \
        --phenoCol $pheno \
        --covarColList $COVAR_LIST \
        --traitType 'binary' \
        --sparseGRMFile $SPARSE_GRM_FILE \
        --sparseGRMSampleIDFile $SPARSE_GRM_ID_FILE \
        --useSparseGRMtoFitNULL TRUE 2>&1 | grep 'Nglmm' | awk -v pheno_var="$pheno" '{print pheno_var "," $2}' >> neff.csv
done

# Process continuous phenotypes
echo "Estimating neff for continuous phenotypes: $cont_phenos"

for pheno in $cont_phenos; do
    echo "Estimating neff for $pheno"
    Rscript extdata/extractNglmm.R \
        --phenoFile $PHENO_FILE \
        --phenoCol $pheno \
        --covarColList $COVAR_LIST \
        --traitType 'quantitative' \
        --sparseGRMFile $SPARSE_GRM_FILE \
        --sparseGRMSampleIDFile $SPARSE_GRM_ID_FILE \
        --useSparseGRMtoFitNULL TRUE 2>&1 | grep 'Nglmm' | awk -v pheno_var="$pheno" '{print pheno_var "," $2}' >> neff.csv
done

echo "Finished estimating neff"

EOF
