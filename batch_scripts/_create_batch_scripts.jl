using Statistics, LinearAlgebra
using Glob
using Printf

number_of_particles = 100_000  # Number of particles to input

particle = "proton"          # "e-" = electrons, "proton" = protons, "gamma" = photons

# Create energy and pitch angle lists
energy_kev_min = 10_000       # Minimum beam energy, keV
energy_kev_max = 1_000_000 #1_000_000       # Maximum beam energy, keV
energy_nbeams = 5             # Number of log-spaced beams to place between minimum and maximum energy
energies_to_simulate = 10.0 .^ LinRange(log10(energy_kev_min), log10(energy_kev_max), energy_nbeams)
energies_to_simulate = round.(energies_to_simulate, digits = 1)

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

  # Don't simulate if we already have data for a given beam
  # TODO

  #=
  if (65 ≤ α ≤ 72) && (E < 200)
    qos = "blanca-lair"
    time_limit = "7-00:00:00"
  end
  =#

  file = open("$(@__DIR__)/$(job_name).sh", "w")
  println(file,
  """
  #!/bin/bash

  #SBATCH --job-name G4EPP_$(job_name)
  #SBATCH --nodes 1
  #SBATCH --ntasks-per-node 40
  #SBATCH --time $(time_limit)
  #SBATCH --output /projects/jucl6426/Aviation_G4EPP/results/log_$(job_name).out
  #SBATCH --qos=$(qos)
  #SBATCH --exclude=bhpc-c5-u7-22,bhpc-c5-u7-23
  #SBATCH --requeue
  #SBATCH --mail-type=ALL
  #SBATCH --mail-user=jucl6426@colorado.edu

  # Terminate on any non-zero exit status
  set -e

  # Run simulation
  cd /projects/jucl6426/Aviation_G4EPP/build/
  ./aviation_G4EPP $(number_of_particles) $(particle) $(energy_string)

  # Copy results to safe folder
  #cp /projects/jucl6426/Aviation_G4EPP/build/results/mlat_65.77deg_input_449.5km_record_450.5km/backscatter_$(input_particle_longname)_input_$(energy_string)keV_deg_$(number_of_particles)particles.csv /projects/jucl6426/G4EPP/results
  #cp /projects/jucl6426/Aviation_G4EPP/build/results/mlat_65.77deg_input_449.5km_record_450.5km/energy_deposition_$(input_particle_longname)_input_$(energy_string)keV_deg_$(number_of_particles)particles.csv /projects/jucl6426/G4EPP/results
  """
  )
  close(file)

  global written += 1
end

println("Wrote $(written) files.")