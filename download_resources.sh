#!/bin/bash

SAIGE_IMAGE=false
PLINK1_9=false
SINGULARITY=false
machine=$(uname)
saige_version="1.1.8"


while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        
        )
        SAIGE_IMAGE=true
        shift
        ;;
        --singularity)
        SINGULARITY=true
        shift
        ;;
        --plink)
        PLINK1_9=true
        shift
        ;;
        *)
        shift
        ;;
    esac
done

if [[ $SAIGE_IMAGE = true ]]; then
    mkdir -p resources/
    if [[ ${SINGULARITY} = true && ! $( test -f "resources/saige-${saige_version}.sif" ) ]]; then
        singularity pull "resources/saige.sif" "docker://wzhou88/saige:${saige_version}"
    elif [[ ${SINGULARITY} = false ]]; then
        docker pull "wzhou88/saige:${saige_version}"
        docker save -o "resources/saige.tar" "wzhou88/saige:${saige_version}"
    fi
fi

if [[ $PLINK1_9 = true ]]; then
    mkdir -p resources/
    if [[ $machine == "Darwin" ]]; then
        echo "Downloading OSX version of plink"
        wget -nc https://s3.amazonaws.com/plink1-assets/plink_mac_20230116.zip --no-check-certificate -P resources/
        unzip -o resources/plink_mac_20230116.zip -d resources/
    elif [[ $machine == "Linux" ]]; then
        echo "Downloading linux version of plink"
        wget -nc https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip --no-check-certificate -P resources/
        unzip -o resources/plink_linux_x86_64_20230116.zip -d resources/
    else
        echo "Operating system not compatible with the code"
    fi
fi
