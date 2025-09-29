#!/bin/bash

#SBATCH --job-name AG4_proton_169047.5keV
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 40
#SBATCH --time 1-00:00:00
#SBATCH --output /projects/jucl6426/Aviation_G4EPP/results/log_AG4_proton_169047.5keV.out
#SBATCH --qos=preemptable
#SBATCH --exclude=bhpc-c5-u7-22,bhpc-c5-u7-23
#SBATCH --requeue
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jucl6426@colorado.edu

# Terminate on any non-zero exit status
set -e

# Run simulation
cd /projects/jucl6426/Aviation_G4EPP/build/
./aviation_G4EPP 100000 proton 169047.5

# Copy results to safe folder
cp /projects/jucl6426/Aviation_G4EPP/build/results/mlat_45deg_input_450km/proton_input_169047.5keV_100000particles_*_spectra.csv /projects/jucl6426/Aviation_G4EPP/results

