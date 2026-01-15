#!/bin/bash

#SBATCH --job-name AG4_alpha_11560505.9keV
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 40
#SBATCH --time 1-00:00:00
#SBATCH --output /projects/jucl6426/Aviation_GLYPHS/results/log_AG4_alpha_11560505.9keV.out
#SBATCH --qos=preemptable
#SBATCH --exclude=bhpc-c5-u7-19,bhpc-c5-u7-22,bhpc-c5-u7-23
#SBATCH --requeue
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jucl6426@colorado.edu

# Terminate on any non-zero exit status
set -e

# Run simulation
cd /projects/jucl6426/Aviation_GLYPHS/build/
./aviation_GLYPHS 100000 alpha 11560505.9

# Copy results to safe folder
cp /projects/jucl6426/Aviation_GLYPHS/build/results/mlat_45deg_input_450km/alpha_input_11560505.9keV_100000particles_*_spectra.csv /projects/jucl6426/Aviation_GLYPHS/results

