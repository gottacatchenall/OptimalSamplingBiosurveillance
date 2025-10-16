using CSV, DataFrames, HTTP, JSON
using CairoMakie
using SpeciesDistributionToolkit
using Statistics
using BiodiversityObservationNetworks

const BONs = BiodiversityObservationNetworks
const SDT = SpeciesDistributionToolkit

include("util.jl")

sdms = read_sdms("SDM_Artifacts_SouthKorea")


uncs = [v[:uncertainty] for (k,v) in sdms]
preds = [v[:prediction] for (k,v) in sdms]

weights = rand(length(uncs))
weights ./= sum(weights)


uncertainty = sum([weights[i]*uncs[i] for i in eachindex(uncs)])
dose = sum([weights[i]*preds[i] for i in eachindex(uncs)])
rescale!(dose)
rescale!(uncertainty)


bon = sample(SimpleRandom(300), uncertainty)

hotspots = sample(AdaptiveHotspot(), bon)

heatmap(uncertainty)
#scatter!([n for n in hotspots.nodes], color=:red)
current_figure()

# Split into quadrants 
p, q = 0.5, 0.5
pred_cutoff = quantile(dose, [p])[begin]
unc_cutoff = quantile(uncertainty, [q])[begin]

A = dose .> pred_cutoff
B = uncertainty .> unc_cutoff


heatmap(nodata(A & B, false), colormap=["#b48ead"], label="Active Surv.")
heatmap!(nodata(.!A & .!B, false), colormap=["#5e81ac"], label="Confirmation")
heatmap!(nodata(.!A & B, false), colormap=["#a3be8c"], label="Discovery")
heatmap!(nodata(A & .!B, false), colormap=["#d08770"], label="Passive Surv.")
axislegend()
current_figure()
