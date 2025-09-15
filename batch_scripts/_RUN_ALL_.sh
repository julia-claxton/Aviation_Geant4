#!/bin/bash

for i in /projects/$USER/G4EPP/source/batch_scripts/*deg.sh; do
  sbatch --quiet "$i"
done
watch -n1 squeue --format=\"%.18i %.15P %.30j %.10u %.10T %.13M %.13l %.8D %R\" --me