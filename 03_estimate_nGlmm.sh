docker run -i -e HOME=${WD} -v ${WD}/:$HOME/ -v /mnt/project/:/mnt/project/ wzhou88/saige:1.3.4 /bin/bash << 'EOF'
for anc in AFR AMR EAS EUR SAS; do 

        sed -i '/setgeno/d' R/SAIGE_extractNeff.R

        # Install the package
        R CMD INSTALL .

        # Define phenotype variables (binary and continuous)
        binary_phenos="Gout Coronary_artery_disease Heart_Failure__HF_ Chronic_obstructive_pulmonary_disease__COPD_ Age_related_macular_degeneration Benign_and_in_situ_intestinal_neoplasms Colon_and_rectum_cancer Inflammatory_bowel_disease Inguinal__femoral__and_abdominal_hernia Interstitial_lung_disease_and_pulmonary_sarcoidosis Non_rheumatic_valvular_heart_disease Pancreatitis Peptic_ulcer_disease Type_2_diabetes Psoriasis Rheumatic_heart_disease Rheumatoid_arthritis Urolithiasis Peripheral_artery_disease Atrial_Fibrillation Varicose_Veins Hypertension Chronic_Renal_Failure Hip_replacement__operation_"

        cont_phenos="BMI Height Alcohol_consumption__drinks_per_week_ Total_cholesterol LDLC HDLC Triglycerides C_reactive_protein__CRP_ Aspartate_aminotransferase__AST_"

        PHENO_FILE="/mnt/project/brava/inputs/phenotypes/brava_with_covariates.tsv"
        COVAR_LIST="age,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
        SPARSE_GRM_FILE="/mnt/project/brava/outputs/step0/brava_${anc}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
        SPARSE_GRM_ID_FILE="/mnt/project/brava/outputs/step0/brava_${anc}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"


        echo "pheno,nglmm" >> Nglmm_$anc

        # Process binary phenotypes
        for pheno in $binary_phenos; do
            Rscript extdata/extractNglmm.R \
                --phenoFile $PHENO_FILE \
                --phenoCol $pheno \
                --covarColList $COVAR_LIST \
                --traitType 'binary' \
                --sparseGRMFile $SPARSE_GRM_FILE \
                --sparseGRMSampleIDFile $SPARSE_GRM_ID_FILE \
                --useSparseGRMtoFitNULL TRUE 2>&1 | grep 'Nglmm' | awk -v pheno_var="$pheno" '{print pheno_var "," $2}' >> Nglmm_$anc
        done

        # Process continuous phenotypes
        for pheno in $cont_phenos; do
            Rscript extdata/extractNglmm.R \
                --phenoFile $PHENO_FILE \
                --phenoCol $pheno \
                --covarColList $COVAR_LIST \
                --traitType 'quantitative' \
                --sparseGRMFile $SPARSE_GRM_FILE \
                --sparseGRMSampleIDFile $SPARSE_GRM_ID_FILE \
                --useSparseGRMtoFitNULL TRUE 2>&1 | grep 'Nglmm' | awk -v pheno_var="$pheno" '{print pheno_var "," $2}' >> Nglmm_$anc
        done

        # Output the results
        mv Nglmm_$anc.csv $HOME/
done
EOF
