using CSV, DataFrames, HTTP, JSON
using CairoMakie
using SpeciesDistributionToolkit
using Statistics
using TernaryDiagrams

const SDT = SpeciesDistributionToolkit

include(joinpath("..", "src", "util.jl"))
include(joinpath("..", "src", "plotting.jl"))

# Load SDMs
sdms = read_sdms(joinpath("artifacts", "India"))

# Get polygons for background land
india_polygon = getpolygon(PolygonData(GADM, Countries), country="IND")
pakistan_polygon = getpolygon(PolygonData(GADM, Countries), country="PAK") 
bng_polygon = getpolygon(PolygonData(GADM, Countries), country="BGD") 
china_polygon = getpolygon(PolygonData(GADM, Countries), country="CHN") 
nepal_polygon = getpolygon(PolygonData(GADM, Countries), country="NPL") 
bhutan_polygon = getpolygon(PolygonData(GADM, Countries), country="BTN") 
afghan_polygon = getpolygon(PolygonData(GADM, Countries), country="AFG") 
myanmar_polygon = getpolygon(PolygonData(GADM, Countries), country="MMR") 
sri_lanka_polygon = getpolygon(PolygonData(GADM, Countries), country="LKA") 
background_land = vcat([pakistan_polygon, bng_polygon, china_polygon, nepal_polygon, bhutan_polygon, afghan_polygon, sri_lanka_polygon, myanmar_polygon]...)

# Get polygons with dashes for contested regions
kashmir_dashes = add_dashes(kashmir_polygon)
ap_dashes = add_dashes(arunuachal_pradesh_polygon)

# Species to use
species = ["Rattus rattus", "Suncus murinus", "Hystrix indica"]

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


# Plotting arguments
background_land_color = :grey92
bbox = SDT.boundingbox(india_polygon)
bbox = bbox.left, bbox.right, bbox.bottom, bbox.top + 2.
disputed_args = (
    color=:grey85, 
    strokewidth=1, 
    strokecolor=:grey50
)
disputed_dashes_args = (
    color=:grey30, 
    strokewidth=1, 
    strokecolor=:grey30
)



# Make Figure
begin 
    fig = Figure(size=(900, 700))

    g = GridLayout(fig[1,1])

    gright = GridLayout(g[1,2])

    ax1 = Axis(
        gright[1,1];
        backgroundcolor=(:red, 0.25),
        ax_args...
    )
    limits!(ax1, bbox...)

    ax2 = Axis(
        gright[2,1];
        backgroundcolor=(:green, 0.25),
        ax_args...
    )
    limits!(ax2, bbox...)


    ax3 = Axis(
        gright[3,1];
        backgroundcolor=(:blue, 0.2),
        ax_args...
    )
    limits!(ax3, bbox...)

    example_weights = [
        [0.8, 0.1, 0.1],
        [0.1, 0.8, 0.1],
        [0.1, 0.1, 0.8]]


    priority1 = quantize(sum([example_weights[1][i]*0.5*(quantize(uncs[i])+quantize(richness[i])) for i in eachindex(uncs)]))
    priority2 = quantize(sum([example_weights[2][i]*0.5*(quantize(uncs[i])+quantize(richness[i])) for i in eachindex(uncs)]))
    priority3 = quantize(sum([example_weights[3][i]*0.5*(quantize(uncs[i])+quantize(richness[i]))  for i in eachindex(uncs)]))


    heatmap!(ax1, priority1, colormap=:lipari)
    heatmap!(ax2, priority2, colormap=:lipari)
    heatmap!(ax3, priority3, colormap=:lipari)



    for ax in [ax1, ax2, ax3]
        poly!(ax, background_land, color=background_land_color)
        poly!(ax, kashmir_polygon; disputed_args...)
        poly!.(ax, kashmir_dashes; disputed_dashes_args...)
        poly!(ax, arunuachal_pradesh_polygon; disputed_args...)
        poly!.(ax, ap_dashes; disputed_dashes_args...)
        lines!(ax, background_land, color=:grey80)
    end 

    ax_ternary = Axis(g[1,1], aspect=1);
    ternaryaxis!(
        ax_ternary; 
        labelx = "",
        labely = "",
        labelz = "",
        tick_fontsize = 14,
    )


    a1 = [w[1] for w in weights]
    a2 = [w[2] for w in weights]
    a3 = [w[3] for w in weights]
    ternaryscatter!(
        ax_ternary, 
        a1,
        a2,
        a3;
        color = [get(Makie.ColorSchemes.magma, w, extrema(maes)) for w in maes],
        marker = :hexagon,
        markersize = 21,
    )

    ew1 = [ew[1] for ew in example_weights]
    ew2 = [ew[2] for ew in example_weights]
    ew3 = [ew[3] for ew in example_weights]

    hidedecorations!(ax_ternary) 
    hidespines!(ax_ternary) 

    scatter!(ax_ternary, [0.15], [0.094], 
        strokewidth=2,
        strokecolor=:white,
        color=:red,
        markersize = 25,
    )
    scatter!(ax_ternary, [0.5], [0.69], 
        strokewidth=2,
        strokecolor=:white,
        color=:blue,
        markersize = 25,
    )
    scatter!(ax_ternary, [0.85], [0.092], 
        strokewidth=2,
        strokecolor=:white,
        color=:green,
        markersize = 25,
    )


    text!(ax_ternary, 0.06, 0.32, text="$(species[1]) weight", rotation=π/3, fontsize=20)
    text!(ax_ternary, 0.32, -0.15, text="$(species[2]) weight", fontsize=20)
    text!(ax_ternary, 0.76, 0.65, text="$(species[3]) weight", rotation=5π/3, fontsize=20)

    colsize!(g, 1, Relative(0.7))

    fig

end 
save("plots/sensitivity.png", fig)
