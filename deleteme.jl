using Statistics, LinearAlgebra
using DelimitedFiles

results_dir = "$(dirname(@__DIR__))/build/results/input_449.5km_record_450.5km"
filenames = readdir(results_dir)
paths = ["$(results_dir)/$(name)" for name in filenames]

modify_times = [stat(file).mtime for file in paths]
_, newest_file_idx = findmax(modify_times)
newest_file = paths[newest_file_idx]

particle = match.(Regex("/input_449.5km_record_450.5km/backscatter_(.*?)_input"), newest_file).captures[1]
n_particles = parse.(Int64, match.(Regex("deg_(.*?)particles.csv"), newest_file).captures[1])
energy = parse.(Int64, match.(Regex("_input_(.*?)keV"), newest_file).captures[1])
pa = parse.(Int64, match.(Regex("keV_(.*?)deg_"), newest_file).captures[1])

backscatter_energy = 0
backscatter_filename = "$(results_dir)/backscatter_$(particle)_input_$(energy)keV_$(pa)deg_$(n_particles)particles.csv"
try
    backscatter_data = readdlm(backscatter_filename, ',', skipstart = 1)
catch
else
    backscatter_data = readdlm(backscatter_filename, ',', skipstart = 1)
    global backscatter_energy = sum(backscatter_data[:,2])
end

deposition_filename = "$(results_dir)/energy_deposition_$(particle)_input_$(energy)keV_$(pa)deg_$(n_particles)particles.csv"
deposition_data = readdlm(deposition_filename, ',', skipstart = 1)
deposition_energy = sum(deposition_data[:,2])

input_energy = energy * n_particles
tracked_energy = backscatter_energy + deposition_energy

println("$(n_particles) $(particle)s @ $(energy/1000) MeV")
println("Tracked $(round(tracked_energy/input_energy, digits = 4)*100)% of input energy")
println("\tBackscatter: $(round(backscatter_energy/input_energy, digits = 4)*100)%")
println("\tDeposition: $(round(deposition_energy/input_energy, digits = 4)*100)%\n")