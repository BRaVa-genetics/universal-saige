run_container () {
  if [[ ${SINGULARITY} = true ]]; then
    singularity exec \
      --env HOME=${WD} \
      --bind ${WD}/:$HOME/ \
      "resources/saige.sif" $cmd
  else
    # Load the Docker image from the tar.gz file
    docker load -i "resources/saige.tar"
    image_id=$(docker images --filter=reference='wzhou88/saige:*' --format "{{.ID}}" | head -n 1)

    docker run \
      -e HOME=${WD} \
      -v ${WD}/:$HOME/ \
      "${image_id}" $cmd
  fi
}