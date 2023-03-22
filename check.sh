#!/bin/bash

function check_container_env() {
  # Check if singularity or docker is installed
  SINGULARITY="$1"
  if [[ ${SINGULARITY} = true ]]; then
    if [[ $(which singularity) == "" ]]; then
      echo "Error: singularity is not installed"
      exit 1
    fi
  else
    if [[ $(which docker) == "" ]]; then
        echo "Error: docker is not installed"
        exit 1
    fi
  fi
  echo "singularity or docker found"
}
