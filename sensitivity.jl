using CSV, DataFrames, HTTP, JSON
using CairoMakie
using SpeciesDistributionToolkit
using Statistics

using TernaryDiagrams

const SDT = SpeciesDistributionToolkit

include("util.jl")

dict = read_sdms("SDM_Artifacts_SouthKorea")

species = ["Rattus norvegicus", "Apodemus peninsulae", "Craseomys regulus"]

uncs = [dict[s][:uncertainty] for s in species]
richness = [dict[s][:prediction] for s in species]

n = 30  # resolution
weights = []
for i in 1:n
    for j in 1:n-i-1
        k = n - i - j 
        push!(weights, (i/n, j/n, k/n))
    end
end

weights


function get_wprime(w, σ)
    w_prime = w .+ σ .* randn(3)
    w_prime .+= minimum(w_prime)
    w_prime ./= sum(w_prime)
    return w_prime
end

σ = 0.01
reps = 100

w = weights[begin]

maes = []
for (wi, w) in enumerate(weights)
    @info "Weight set $wi/$(length(weights))"
    priority = sum([w[i]*uncs[i] for i in eachindex(w)])

    mae_sum = 0.
    for _ in 1:reps
        w_prime = get_wprime(w, σ)
        priority_prime = sum([w_prime[i]*uncs[i] for i in eachindex(w)])
        mae = mean(abs.(priority - priority_prime))
        mae_sum += mae
    end 
    avg_mae = mae_sum / reps
    push!(maes, avg_mae)
end

maes



fig = Figure(size=(1000,1000))

ax1 = Axis(fig[1,1], leftspinecolor=:red, rightspinecolor=:red, topspinecolor=:red, bottomspinecolor=:red,  spinewidth=6)
ax2 = Axis(fig[1,2], leftspinecolor=:green, rightspinecolor=:green, topspinecolor=:green, bottomspinecolor=:green,  spinewidth=6)
ax3 = Axis(fig[2,1], leftspinecolor=:blue, rightspinecolor=:blue, topspinecolor=:blue, bottomspinecolor=:blue,  spinewidth=6)


example_weights = [
    [0.8, 0.1, 0.1],
    [0.1, 0.8, 0.1],
    [0.1, 0.1, 0.8]]


priority1 = quantize(sum([example_weights[1][i]*0.5*(quantize(uncs[i])+quantize(richness[i])) for i in eachindex(uncs)]))
priority2 = quantize(sum([example_weights[2][i]*0.5*(quantize(uncs[i])+quantize(richness[i])) for i in eachindex(uncs)]))
priority3 = quantize(sum([example_weights[3][i]*0.5*(quantize(uncs[i])+quantize(richness[i]))  for i in eachindex(uncs)]))


heatmap!(ax1, priority1)
heatmap!(ax2, priority2)
heatmap!(ax3, priority3)


ax_ternary = Axis(fig[2,2]);
ternaryaxis!(
    ax_ternary; 
    labelx = "w1",
    labely = "w2",
    labelz = "w3",

)
a1 = [w[1] for w in weights]
a2 = [w[2] for w in weights]
a3 = [w[3] for w in weights]
ternaryscatter!(
    ax_ternary, 
    a1,
    a2,
    a3;
    color = [get(Makie.ColorSchemes.Spectral, w, extrema(maes)) for w in maes],
    marker = :hexagon,
    markersize = 25,
)

ew1 = [ew[1] for ew in example_weights]
ew2 = [ew[2] for ew in example_weights]
ew3 = [ew[3] for ew in example_weights]
ternaryscatter!(
    ax_ternary, 
    ew1,
    ew2,
    ew3;
    color = [:red, :green, :blue],
    marker = :star,
    markersize = 25,
)
hidedecorations!(ax_ternary) 
hidespines!(ax_ternary) 

fig

save("sensitivity.png", fig)



