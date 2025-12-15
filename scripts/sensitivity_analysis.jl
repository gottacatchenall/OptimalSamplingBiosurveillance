@info "Running Sensitivity Analysis"

using CSV, DataFrames, JSON
using CairoMakie
using SpeciesDistributionToolkit
using Statistics

const AG = SpeciesDistributionToolkit.SimpleSDMPolygons.AG
const SDT = SpeciesDistributionToolkit

include(joinpath("..", "src", "util.jl"))
include(joinpath("..", "src", "plotting.jl"))

# Load SDMs
sdms = read_sdms(joinpath("artifacts", "India"))

# Species to use
species = ["Rattus rattus", "Suncus murinus", "Hystrix indica"]

# Get uncertainties for each species
uncs = [sdms[s][:uncertainty] for s in species]
richness = [sdms[s][:prediction] for s in species]

# Create weights that span all possible values
n = 40  # Resolution of ternary plot points
weights = []
for i in 1:n
    for j in 1:n-i-1
        k = n - i - j 
        push!(weights, (i/n, j/n, k/n))
    end
end

# Renormalize nudged weights
function get_wprime(w, σ)
    w_prime = w .+ σ .* randn(3)
    w_prime .+= minimum(w_prime)
    w_prime ./= sum(w_prime)
    return w_prime
end

# Noise amount
σ = 0.01
# Replicates for unique weight set
reps = 100

# Compute average difference between priority computed with weights before and after nudging them with noise
maes = []
for (wi, w) in enumerate(weights)
    @info "\tWeight set $wi/$(length(weights))"
    priority = sum([w[i]*(uncs[i] + richness[i]) for i in eachindex(w)])

    mae_sum = 0.
    for _ in 1:reps
        w_prime = get_wprime(w, σ)
        priority_prime = sum([w_prime[i]*(uncs[i]+richness[i]) for i in eachindex(w)])
        mae = mean(abs.(priority - priority_prime))
        mae_sum += mae
    end 
    avg_mae = mae_sum / reps
    push!(maes, avg_mae)
end


# Write Results
df = DataFrame(mae=[], rattus_weight=[], suncus_weight=[], hystrix_weight = [])
for i in eachindex(maes)
    push!(df, (maes[i], weights[i]...))
end 
mkpath("artifacts")
CSV.write(joinpath("artifacts", "sensitivity.csv"), df)
