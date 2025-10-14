using Statistics, LinearAlgebra
using Glob
using Printf

include("/Users/luna/Research/Aviation_Radiation/code/SpectrumParser.jl")
results_dir = "/Users/luna/Research/geant4/Aviation_GLYPHS/results"
beam_particles, beam_energies_keV = get_beams(results_dir)

number_of_particles = 100_000  # Number of particles to input

particle = "e-" #"proton"          # "e-" = electrons, "proton" = protons, "gamma" = photons

# Create energy and pitch angle lists
energy_kev_min = 100_000       # Minimum beam energy, keV
energy_kev_max = 100_000_000   # Maximum beam energy, keV
energy_nbeams = 150            # Number of log-spaced beams to place between minimum and maximum energy
energies_to_simulate = 10.0 .^ LinRange(log10(energy_kev_min), log10(energy_kev_max), energy_nbeams)
energies_to_simulate = round.(energies_to_simulate, digits = 1)










@warn "undo e- changes"
energies_to_simulate = 10_000 # XXX TEMP











# Create shell scripts
rm.(glob("*keV.sh", @__DIR__))
written = 0
skipped = 0

for E in energies_to_simulate
  input_particle_longname = particle == "e-" ? "electron" : particle
  energy_string = @sprintf "%.1f" E
  job_name = "AG4_$(input_particle_longname)_$(energy_string)keV"
  qos = "preemptable"
  time_limit = "1-00:00:00"

  # Send long-runtime beams to blanca-lair
  if E ≥ 50_000_000
    qos = "blanca-lair"
    time_limit = "7-00:00:00"
  end
  
  #=
  # Don't simulate if we already have data for a given beam
  if E in beam_energies_keV
    global skipped += 1
    continue
  end
  =#

  file = open("$(@__DIR__)/$(job_name).sh", "w")
  println(file,
  """
  #!/bin/bash

  #SBATCH --job-name $(job_name)
  #SBATCH --nodes 1
  #SBATCH --ntasks-per-node 40
  #SBATCH --time $(time_limit)
  #SBATCH --output /projects/jucl6426/Aviation_GLYPHS/results/log_$(job_name).out
  #SBATCH --qos=$(qos)
  #SBATCH --exclude=bhpc-c5-u7-19,bhpc-c5-u7-22,bhpc-c5-u7-23
  #SBATCH --requeue
  #SBATCH --mail-type=ALL
  #SBATCH --mail-user=jucl6426@colorado.edu

  # Terminate on any non-zero exit status
  set -e

  # Run simulation
  cd /projects/jucl6426/Aviation_GLYPHS/build/
  ./aviation_GLYPHS $(number_of_particles) $(particle) $(energy_string)

  # Copy results to safe folder
  cp /projects/jucl6426/Aviation_GLYPHS/build/results/mlat_45deg_input_450km/$(input_particle_longname)_input_$(energy_string)keV_$(number_of_particles)particles_*_spectra.csv /projects/jucl6426/Aviation_GLYPHS/results
  """
  )
  close(file)

  global written += 1
end

println("Wrote $(written) files, skipped $(skipped) files.")