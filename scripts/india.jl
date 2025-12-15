@info "Running India case study..."


using ColorBlendModes
using Random

include(joinpath("..", "src", "sdms.jl"))
include(joinpath("..", "src", "plotting.jl"))
include(joinpath("..", "src", "util.jl"))
using BiodiversityObservationNetworks
const AG = SDT.SimpleSDMPolygons.AG


# Load historical sampling sites from Rodrigues et al. (1978)
historical_coords = unique([(r.longitude, r.latitude) for r in eachrow(CSV.read(joinpath("data", "lassa_india.csv"), DataFrame))])

# Download occurrences from GBIF and group by species
doi = "10.15468/dl.b3z33h"
india_occurrences = GBIF.download(doi)
occurrence_by_species = group_occurrences_by_species(india_occurrences)

india_polygon = getpolygon(PolygonData(GADM, Countries), country="IND") 

# Remove contested areas from India polygon
kashmir_polygon = india_polygon[2]
arunuachal_pradesh_polygon = india_polygon[5]
india_polygon = india_polygon - kashmir_polygon - arunuachal_pradesh_polygon

# Create dashed polygons for contested areas
kashmir_dashes = add_dashes(kashmir_polygon)
ap_dashes = add_dashes(arunuachal_pradesh_polygon)

# Get polygon for background land
pakistan_polygon = getpolygon(PolygonData(GADM, Countries), country="PAK") 
bng_polygon = getpolygon(PolygonData(GADM, Countries), country="BGD") 
china_polygon = getpolygon(PolygonData(GADM, Countries), country="CHN") 
nepal_polygon = getpolygon(PolygonData(GADM, Countries), country="NPL") 
bhutan_polygon = getpolygon(PolygonData(GADM, Countries), country="BTN") 
afghan_polygon = getpolygon(PolygonData(GADM, Countries), country="AFG") 
myanmar_polygon = getpolygon(PolygonData(GADM, Countries), country="MMR") 
sri_lanka_polygon = getpolygon(PolygonData(GADM, Countries), country="LKA") 
background_land = vcat([pakistan_polygon, bng_polygon, china_polygon, nepal_polygon, bhutan_polygon, afghan_polygon, sri_lanka_polygon, myanmar_polygon]...)
background_land_color = :grey92



# Get Environmental Predictors
environmental_layers = [Float32.(SDMLayer(RasterData(WorldClim2, BioClim); resolution=2.5, layer=i, SDT.boundingbox(india_polygon)...)) for i in 1:19]
mask!(environmental_layers, india_polygon)

# Fit SDMs
@info "\tFitting SDMs for India Case Study..."
sdms = Dict([s=>fit_sdm(occurrence_by_species[s], environmental_layers) for s in keys(occurrence_by_species)])
write_sdm_artifacts(joinpath("artifacts", "India"), sdms)

# Load SDMs (you can skip previous 2 lines if they already exist)
sdms = read_sdms(joinpath("artifacts", "India"))

weights = [0.11, 0.34, 0.2, 0.35]

# Even Dose vs. Uncertainty Weighting
w_dose = 0.5
w_uncertainty = 1 - w_dose

# Create weighted dose and uncertainty components of priority
weighted_uncertainty = sum([weights[i]*quantize(r[:uncertainty]) for (i,r) in enumerate(values(sdms))])
weighted_dose = sum([weights[i]*quantize(r[:prediction]) for (i,r) in enumerate(values(sdms))])

# Create overall priority map
priority = rescale(w_dose .* weighted_dose + w_uncertainty .* weighted_uncertainty, (0,1))

function get_num_species_in_prevalence(sdms)
    species_in_prevalence = similar(first(sdms)[2][:prediction])
    species_in_prevalence.grid .= 0

    for (s, sdm) in sdms
        P = sdm[:prediction]
        Q = Bool.(discretize(quantize(P, 2), 2) .- 1)
        species_in_prevalence.grid[findall(Q)] .+= 1
    end
    species_in_prevalence
end

function get_num_species_in_discovery(sdms)
    species_in_discovery = similar(first(sdms)[2][:prediction])
    species_in_discovery.grid .= 0
    for (s, sdm) in sdms
        U = sdm[:uncertainty]
        Q = Bool.(discretize(quantize(U, 2), 2) .- 1)
        species_in_discovery.grid[findall(Q)] .+= 1
    end
    species_in_discovery
end

function get_most_important_host()
    hosts = collect(keys(sdms))
    most_important_host = similar(priority)

    for i in eachindex(most_important_host)
        scores = [
            weights[h] * 
                (w_dose * sdms[hosts[h]][:prediction][i] + 
                w_uncertainty * sdms[hosts[h]][:uncertainty][i])
            for h in eachindex(hosts)
        ]
        most_important_host.grid[i] = findmax(scores)[2]
    end
    return hosts, most_important_host
end



# Create dose vs. uncertainty bivariate map
nbreaks1 = 4
bivar, colormatrix = get_bivariate(weighted_dose, weighted_uncertainty; nbreaks=nbreaks1)


# Create bivariate map of number of species in both prevalence and discovery sampling
prev, disc = get_num_species_in_prevalence(sdms), get_num_species_in_discovery(sdms)
nbreaks2 = 4
bivar2, colormatrix2 = get_bivariate(
    prev, 
    disc; 
    nbreaks=nbreaks2,
    high1=colorant"#b269a9", 
    high2=colorant"#78c6c7",
)

# Sample a BON
Random.seed!(42)
bon = BiodiversityObservationNetworks.sample(
    BalancedAcceptance(), 
    priority, 
    inclusion=BiodiversityObservationNetworks.tilt(priority, 5)
)


# Get most important host at each site
hosts, most_important_host = get_most_important_host()


# Arguments for plotting
ax_settings = (
    ;
    aspect = 1,
    xgridvisible=false,
    ygridvisible=false
)

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

most_impt_colors = [colorant"#a3be8c", colorant"#b48ead", colorant"#88c0d0", colorant"#5e81ac"]

CairoMakie.activate!(; px_per_unit=3)

# Get padded boundingbox for axis limits
bbox = SDT.boundingbox(india_polygon)
bbox = bbox.left, bbox.right, bbox.bottom, bbox.top + 2.


@info "\tMaking Figure for India Case Study..."
begin 
    fig = Figure(size=(900,900))
    ax_dose = Axis(fig[1,1]; ax_settings...)
    limits!(ax_dose, bbox...)

    ax_uncertainty = Axis(fig[1,2]; aspect=1, ax_settings...)
    limits!(ax_uncertainty, bbox...)

    ax_priority = Axis(fig[2,1]; ax_settings...)
    limits!(ax_priority, bbox...)

    ax_most_impt_host = Axis(fig[2,2]; ax_settings...)
    limits!(ax_most_impt_host, bbox...)

    poly!(ax_dose, background_land, color=background_land_color)
    lines!(ax_dose, background_land, color=:grey80)
    heatmap!(ax_dose, bivar, colormap=vec(colormatrix))
    poly!(ax_dose, kashmir_polygon; disputed_args...)
    poly!.(ax_dose, kashmir_dashes; disputed_dashes_args...)
    poly!(ax_dose, arunuachal_pradesh_polygon; disputed_args...)
    poly!.(ax_dose, ap_dashes; disputed_dashes_args...)
    text!(ax_dose, 94, 32, text="A", fontsize=30, font=:bold)

    ax_inset = Axis(fig[1, 1],
        width=Relative(0.28),
        height=Relative(0.28),
        aspect=1,
        halign=0.82,
        valign=0.1,
        xticksvisible=false,
        yticksvisible=false,
        xticklabelsvisible=false,
        yticklabelsvisible=false,
        xlabel = "",
        ylabel = "",
    )
    heatmap!(ax_inset, (reshape(colormatrix, (nbreaks1,nbreaks1))))

    annotation!(ax_dose, -83, 0, 93, 8; text = "Dose", style = Ann.Styles.LineArrow(head=Ann.Arrows.Head(length=9)))
    annotation!(ax_dose, 0, -55, 84.7, 17.5; text = " ", style = Ann.Styles.LineArrow(head=Ann.Arrows.Head(length=9)))
    text!(
        ax_dose,
        [(85.1, 8.9)],
        text = "Uncertainty",
        rotation = π/2,
    )
    
    poly!(ax_uncertainty, background_land, color=background_land_color)
    lines!(ax_uncertainty, background_land, color=:grey80)
    heatmap!(ax_uncertainty, bivar2, colormap=vec(colormatrix2))
    poly!(ax_uncertainty, kashmir_polygon; disputed_args...)
    poly!.(ax_uncertainty, kashmir_dashes; disputed_dashes_args...)
    poly!(ax_uncertainty, arunuachal_pradesh_polygon; disputed_args...)
    poly!.(ax_uncertainty, ap_dashes; disputed_dashes_args...)
    text!(ax_uncertainty, 94, 32, text="B", fontsize=30, font=:bold)
    ax_inset = Axis(fig[1, 2],
        width=Relative(0.25),
        height=Relative(0.25),
        aspect=1,
        halign=0.81,
        valign=0.18,
        xticks=1:4,
        yticks=1:4,
        xticklabelsize=10,
        yticklabelsize=10,
        xticksvisible=false,
        yticksvisible=false,
        xlabel = "# Prevalence Sp.",
        ylabel = "# Discovery Sp."
    )
    heatmap!(ax_inset, (reshape(colormatrix2, (nbreaks2,nbreaks2))))

    poly!(ax_priority, background_land, color=background_land_color)
    lines!(ax_priority, background_land, color=:grey80)
    heatmap!(ax_priority, priority, colormap=:lipari)
    scatter!(ax_priority, [n for n in bon], color=colorant"#fff", strokecolor=:black, strokewidth=1)
    scatter!(ax_priority, historical_coords, color=colorant"#9effff", marker=:x, markersize=12, strokecolor=:black, strokewidth=1)
    poly!(ax_priority, kashmir_polygon; disputed_args...)
    poly!.(ax_priority, kashmir_dashes; disputed_dashes_args...)
    poly!(ax_priority, arunuachal_pradesh_polygon; disputed_args...)
    poly!.(ax_priority, ap_dashes; disputed_dashes_args...)
   
    cbar_axis = Axis(
        fig[2,1],
        width=Relative(0.4),
        height=Relative(0.03),
        halign=0.8,
        valign=0.09,
        xticksvisible=false,
        yticksvisible=false,
        xticklabelsvisible=false,
        yticklabelsvisible=false,
        title="Sampling Priority",
        titlesize = 12,
        titlefont=:regular,
        titlealign=:left
    )
    cmap = [get(Makie.ColorSchemes.lipari, i) for i in 0:0.025:1]
    X = hcat(0:0.025:1)
    heatmap!(cbar_axis, X, colormap=cmap)
    annotation!(ax_priority, -45, 0, 92.3, 10.9; text = " ", style = Ann.Styles.LineArrow())



    text!(ax_priority, 94, 32, text="C", fontsize=30, font=:bold)

    poly!(ax_most_impt_host, background_land, color=background_land_color)
    lines!(ax_most_impt_host, background_land, color=:grey80)
    heatmap!(ax_most_impt_host, most_important_host, colormap=most_impt_colors)
    poly!(ax_most_impt_host, kashmir_polygon; disputed_args...)
    poly!.(ax_most_impt_host, kashmir_dashes; disputed_dashes_args...)
    poly!(ax_most_impt_host, arunuachal_pradesh_polygon; disputed_args...)
    poly!.(ax_most_impt_host, ap_dashes; disputed_dashes_args...)
    text!(ax_most_impt_host, 94, 32, text="D", fontsize=30, font=:bold)


    axislegend(
        ax_most_impt_host,
        [MarkerElement(color = c, marker = :rect, markersize=23) for c in most_impt_colors],
        hosts,
        "Largest Host Contributor",
        position=:rb
    )

    fig
end

save(joinpath("plots", "india.png"), fig)


