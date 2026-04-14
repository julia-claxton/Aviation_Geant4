using Statistics, LinearAlgebra
using Glob
using Printf

include("/Users/luna/Research/Aviation_Radiation/code/SpectrumParser.jl")
results_dir = "/Users/luna/Research/Aviation_Radiation/data/GLYPHS"
existing_particles, existing_energies = get_beams(results_dir)

number_of_particles = 100_000  # Number of particles to input
alpha_energies = 1e6 .* [1.467029285016976, 1.7817875218749928, 2.1434580789902244, 2.5531026851640726, 3.0195499666152683, 3.5436236292778935, 4.125627948346392, 4.774463797977579, 5.490273042738451, 6.282166640382038, 7.159519110086783, 8.122386626142525, 9.170597549523674, 10.313477444726429, 11.560505948882511, 12.911525353931168, 14.376071249205166, 15.963767866189537, 17.674457592576562, 20.53486106783365, 26.080040217819565, 35.821267113820014, 52.667749722531134, 78.32854513709988, 114.50220914157535, 165.98409421382254]
proton_energies = 1e6 .* [0.4919968561663142, 0.6202564919316679, 0.7632187244759833, 0.9246328837346689, 1.1043104921012188, 1.3019751267372308, 1.521935705828299, 1.7640765698686003, 2.032968417343461, 2.3285487124854978, 2.650701655257517, 3.0041418809864373, 3.388794404104349, 3.809461740784525, 4.271004102390247, 4.773387063190187, 5.316554218799155, 5.9053912772825905, 6.544811688653951, 7.234780943620419, 7.98022843236318, 8.78609801629206, 9.652363544291253, 11.097336568805867, 13.890406704078622, 18.78300886683073, 27.226294269353872, 40.07138487336281, 58.168090056375284, 83.9158247545856]
electron_energies = [63.245540618896484, 97.97958374023438, 138.5640869140625, 183.30308532714844, 238.11758422851562, 305.20489501953125, 385.16229248046875, 520.48046875, 752.9939575195312, 1081.665283203125, 1529.7060546875, 2121.3203125, 2893.960205078125, 3728.6064453125, 4906.12060546875, 6500.0]

alpha_energies = round.(alpha_energies, digits = 1)
proton_energies = round.(proton_energies, digits = 1)
electron_energies = round.(electron_energies, digits = 1)

energies = [proton_energies..., alpha_energies..., electron_energies...]
particles = [repeat(["proton"], length(proton_energies))..., repeat(["alpha"], length(alpha_energies))..., repeat(["e-"], length(electron_energies))...]

# Create shell scripts
rm.(glob("*keV.sh", @__DIR__))
written = 0
skipped = 0
lair = 0

for i in eachindex(energies)
  particle = particles[i]
  E = energies[i]
  input_particle_longname = particle == "e-" ? "electron" : particle
  energy_string = @sprintf "%.1f" E
  job_name = "AG4_$(input_particle_longname)_$(energy_string)keV"
  qos = "preemptable"
  time_limit = "1-00:00:00"

  # Don't simulate if we already have data for a given beam
  if E in existing_energies
    global skipped += 1
    continue
  end

  # Send long-runtime beams to blanca-lair
  if particle == "alpha"
    qos = "blanca-lair"
    time_limit = "3-00:00:00"
    global lair += 1
  end

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
  #SBATCH --mail-user=julia.claxton@colorado.edu

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

println("$(written) files written.")
println("$(skipped) files skipped.")
println("$(lair) files sent to blanca-lair.")