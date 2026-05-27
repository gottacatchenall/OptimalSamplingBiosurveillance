@info "Running South Korea case study..."

using BiodiversityObservationNetworks
using Random
using JSON
const BONs = BiodiversityObservationNetworks

Random.seed!(123)

# Load Dependencies 
include(joinpath("..", "src", "sdms.jl"))
include(joinpath("..", "src", "plotting.jl"))
include(joinpath("..", "src", "util.jl"))
include(joinpath("..", "src", "kriging.jl"))

# Hosts and Viruses of interest
hosts = ["Rattus norvegicus", "Apodemus agrarius", "Crocidura lasiura"]
viruses = ["Orthohantavirus hantanense", "Orthohantavirus seoulense"]


# Download Occurrences from GBIF
doi = "10.15468/dl.wzcxp9"
occurrence_records = GBIF.download(doi)

# Filter occurrences without coordinate (or with a coordinate in the Atlantic)
filter!(x -> !ismissing(place(x)), elements(occurrence_records))
filter!(x -> longitudes(x) > 100, elements(occurrence_records))
occurrence_by_species = Dict([s=>Occurrences(filter(x -> startswith(s)(entity(x)), elements(occurrence_records))) for s in hosts])

# Get Polygons
south_korea_polygon = getpolygon(PolygonData(GADM, Countries), country="KOR")
north_korea_polygon = getpolygon(PolygonData(GADM, Countries), country="PRK")

# Load and mask CHELSA Bioclim layers
environmental_layers = [Float32.(SDMLayer(RasterData(CHELSA2, BioClim); layer=i, SDT.boundingbox(south_korea_polygon)...)) for i in 1:19]
mask!(environmental_layers, south_korea_polygon)

# Fit and Write SDMs
@info "\tFitting SDMs for South Korea Case Study..."
sdms = Dict([
    s=> fit_sdm(
        occurrence_by_species[s], 
        environmental_layers, 
        pseudoabsence_buffer_distance=8.
    ) for s in hosts]
)
write_sdm_artifacts(joinpath("artifacts", "SouthKorea"), sdms)

# Load SDMs (if already fit, previous 2 lines can be skipped)
sdms = read_sdms(joinpath("artifacts", "SouthKorea"))

# Number of breaks in bivariate plots
nbreaks = 3

# Bounding box for plots
bounds = (125.8, 129.8, 33., 38.6)



# ------------------------------------------------------------------------------------
# Kriging for Prevalence 
# ------------------------------------------------------------------------------------
@info "\tRunning kriging..."
hv_pairs = get_prevalence_data()
interpolated_prev = []
for p in hv_pairs
    gp_result, gp_var = fit_gp(p[:coordinates], p[:prevalence], environmental_layers[begin], south_korea_polygon)
    push!(
        interpolated_prev,
        Dict(
            :host=>p[:host],
            :virus=>p[:virus],
            :prev_predict=>gp_result,
            :prev_var=>gp_var,
        )
    )
end



# ------------------------------------------------------------------------------------
# Set Weights
# ------------------------------------------------------------------------------------
Random.seed!(42)
host_weights = rand(Dirichlet([2. for _ in 1:length(hosts)]))
virus_weights = rand(Dirichlet([1. for _ in 1:length(viruses)]))
w_prevalence = 0.5
w_dose = 0.5
w_uncertainty = 1 - w_dose

dose_priority = similar(first(sdms)[2][:prediction])
dose_priority.grid .= 0

uncertainty_priority = similar(first(sdms)[2][:prediction])
uncertainty_priority.grid .= 0

prevalence_priority = similar(first(environmental_layers))
prevalence_priority.grid .= 0


# ------------------------------------------------------------------------------------
# Get weights for each (host, virus) pair with prevalence data
# ------------------------------------------------------------------------------------

rescaled_hv_weights = zeros(length(hosts), length(viruses))
for hi in eachindex(hosts), vi in eachindex(viruses)
    idx = findfirst([p[:host] == hosts[hi] && p[:virus] == viruses[vi] for p in interpolated_prev])
    if !isnothing(idx)
        rescaled_hv_weights[hi,vi] = host_weights[hi] * virus_weights[vi]
    end
end
rescaled_hv_weights ./= sum(rescaled_hv_weights)

# ------------------------------------------------------------------------------------
# Compute Prevalence Priority 
# ------------------------------------------------------------------------------------
for hi in eachindex(hosts)
    w_host = host_weights[hi]
    P = quantize(sdms[hosts[hi]][:prediction])
    U = quantize(sdms[hosts[hi]][:uncertainty])
    global dose_priority += w_host .* (w_dose .* P)   
    global uncertainty_priority += w_host * (w_uncertainty .* U)
    for vi in eachindex(viruses)
        # Does this have a spatial prev estimate?
        idx = findfirst([p[:host] == hosts[hi] && p[:virus] == viruses[vi] for p in interpolated_prev])
        if !isnothing(idx)
            prev_unc = quantize(interpolated_prev[idx][:prev_var])
            this_pair = rescaled_hv_weights[hi, vi] .* prev_unc
            global prevalence_priority += this_pair
        end 
    end
end

# Idk why the coordinates get offset by like 10^-10 
prevalence_priority.x = dose_priority.x
prevalence_priority.y = dose_priority.y

# ------------------------------------------------------------------------------------
# Stratification:  
# ------------------------------------------------------------------------------------
# Top 50% of Dose and Top 50% of Prevalence Uncertainty: Prevalence Sampling
# Top 50% of Uncertainty: Occurrence Sampling

top_dose = Bool.(discretize(quantize(dose_priority), 2)  .- 1)
top_prev = Bool.(discretize(quantize(prevalence_priority, 2), 2) .- 1)
top_unc = Bool.(discretize(quantize(uncertainty_priority, 2), 2) .- 1)

prev_sampling_region = top_dose .& top_prev
occ_sampling_region = top_unc 
both_sampling_region = prev_sampling_region .& occ_sampling_region

OCCURRENCE_SAMPLING = 1
PREVALENCE_SAMPLING = 2
BOTH_SAMPLING = 3

strata = Int64.(similar(top_dose))
strata.grid .= 0
strata.grid[findall(top_unc)] .= OCCURRENCE_SAMPLING
strata.grid[findall(prev_sampling_region)] .= PREVALENCE_SAMPLING 
strata.grid[findall(prev_sampling_region .& occ_sampling_region)] .= BOTH_SAMPLING


# Compute priority within occurrence priority regions
within_occ_priority = Float64.(deepcopy(occ_sampling_region))
nodata!(within_occ_priority, 0)
within_occ_priority.grid .= (dose_priority + uncertainty_priority).grid
quantize!(within_occ_priority)

# Compute priority within prevalence priority regions
within_prev_priority = Float64.(deepcopy(prev_sampling_region))
nodata!(within_prev_priority, 0)
within_prev_priority.grid .= (dose_priority + prevalence_priority).grid
quantize!(within_prev_priority)



# ------------------------------------------------------------------------------------
# Get total priority 
# ------------------------------------------------------------------------------------
function get_total_priority(w_prevalence, within_prev_priority, within_occ_priority)
    total_priority = similar(within_prev_priority)
    total_priority.grid .= 0
    total_priority.indices .= within_prev_priority.indices .| within_occ_priority.indices
    total_priority.grid[total_priority.indices] += (w_prevalence .* within_prev_priority.grid[total_priority.indices])
    total_priority.grid[total_priority.indices] += ((1-w_prevalence) .* within_occ_priority.grid[total_priority.indices])
    quantize!(total_priority)
    return total_priority
end 

w_prevalence = 0.5

total_priority = get_total_priority(w_prevalence, within_prev_priority, within_occ_priority)

tilt(x, α) = exp.(α .* x)

Random.seed!(420)
bon = BiodiversityObservationNetworks.sample(
    BalancedAcceptance(), 
    total_priority, 
    inclusion=tilt(rescale(total_priority, (0,1)), 2)
)


# ------------------------------------------------------------------------------------
# Make Main South Korea Figure  
# ------------------------------------------------------------------------------------

markers = [:utriangle, :circle, :x,]
markersizes = [18, 15, 20]

marker = [markers[strata[n]]  for n in bon]
markersize = [markersizes[strata[n]]  for n in bon]
arrowcolor = :grey40


@info "\tVisualizing South Korea Case Study..."
begin 

f = Figure(size=(1100, 700))
g_base = GridLayout(f[1,1])

g_top_left = GridLayout(g_base[1,1])
g_top_right = GridLayout(g_base[1,2])
#g_bottom = GridLayout(g_base[2,:])

# -------------- Dose vs. Dose Uncertainty Bivariate  --------------
bivar, colormatrix = get_bivariate(dose_priority, uncertainty_priority; nbreaks=nbreaks)
ax = Axis(g_top_left[1,1], aspect=DataAspect(), xgridvisible=false, ygridvisible=false)
limits!(ax, bounds...)
text!(ax, 129.2, 37.95, text="A", fontsize=25, font=:bold, color=:grey20)
hidedecorations!(ax)
heatmap!(ax, bivar, colormap=vec(colormatrix))
poly!(ax, north_korea_polygon, color=:grey85, strokewidth=2, strokecolor=:grey40)
add_bivariate_legend!(g_top_left[1,1], colormatrix, nbreaks, xlabel="", ylabel="", halign=0.86, width=0.26)
annotation!(ax, -45, 0, 129.7, 33.15; text = "Dose", fontsize=11, style = Ann.Styles.LineArrow(head = Ann.Arrows.Head(length=7)))
annotation!(ax, 0, -70, 128.3, 34.5; text = " ", fontsize=11, style = Ann.Styles.LineArrow(head = Ann.Arrows.Head(length=7)))
text!(ax, 128.13, 33.3, text="Dose\nUncertainty", fontsize=11, rotation=π/2)


# -------------- Dose vs. Prevalence Uncertainty Bivariate  --------------
bivar, colormatrix = get_bivariate(dose_priority, prevalence_priority; nbreaks=nbreaks,  high2=colorant"#f2aa7a")
ax = Axis(g_top_left[1,2], aspect=DataAspect(), xgridvisible=false, ygridvisible=false)
limits!(ax, bounds...)
text!(ax, 129.2, 37.95, text="B", fontsize=25, font=:bold, color=:grey20)
hidedecorations!(ax)
heatmap!(ax, bivar, colormap=vec(colormatrix))
poly!(ax, north_korea_polygon, color=:grey85, strokewidth=2, strokecolor=:grey40)
add_bivariate_legend!(g_top_left[1,2], colormatrix, nbreaks, xlabel="", ylabel="", halign=0.86, width=0.26)
annotation!(ax, -45, 0, 129.7, 33.15; text = "Dose", fontsize=11, style = Ann.Styles.LineArrow(head = Ann.Arrows.Head(length=7)))
annotation!(ax, 0, -70, 128.3, 34.5; text = " ", fontsize=11, style = Ann.Styles.LineArrow(head = Ann.Arrows.Head(length=7)))
text!(ax, 128.13, 33.3, text="Prevalence\nUncertainty", fontsize=11, rotation=π/2)


# -------------- Sampling Type Strata  --------------
stratum_cmap = [ colorant"#5e81ac", colorant"#bf616a",colorant"#ebcb8b"]

ax = Axis(g_top_left[2,1], aspect=DataAspect(), xgridvisible=false, ygridvisible=false)
limits!(ax, bounds...)
text!(ax, 129.2, 37.95, text="C", fontsize=25, font=:bold, color=:grey20)
hidedecorations!(ax)
heatmap!(ax, nodata(strata, !iszero), colormap=[:grey90])
poly!(ax, north_korea_polygon, color=:grey90, strokewidth=2, strokecolor=:grey40)
heatmap!(ax, nodata(strata, 0), colormap=stratum_cmap)
leg = axislegend(
    ax,
    [MarkerElement(color = c, marker = :rect, markersize=23) for c in stratum_cmap],
    ["Occurrence", "Prevalence", "Both"],
    position=:rb,
    labelsize=12,
    rowgap=1,
    framevisible=false
)
f 


# -------------- Discovery vs. Prevalence Priority  --------------
ax = Axis(g_top_left[2,2], aspect=DataAspect(), xgridvisible=false, ygridvisible=false)
limits!(ax, bounds...)
text!(ax, 129.2, 37.95, text="D", fontsize=25, font=:bold, color=:grey20)
hidedecorations!(ax)
heatmap!(ax, environmental_layers[1], colormap=[:grey90])
heatmap!(ax, quantize(within_prev_priority), colormap=:Purples)
heatmap!(ax, quantize(within_occ_priority), colormap=:Blues)
poly!(ax, north_korea_polygon, color=:grey85, strokewidth=2, strokecolor=:grey40)
add_colorbar!(g_top_left[2,2], Makie.ColorSchemes.Purples, halign=0.85, valign=0.05, width=0.4, title="Prevalence Priority")
add_colorbar!(g_top_left[2,2], Makie.ColorSchemes.Blues, halign=0.85, width=0.4, title="Discovery Priority")

# -------------- Total Priority  --------------

ax = Axis(g_top_right[1,1], aspect=DataAspect(), xgridvisible=false, ygridvisible=false)
limits!(ax, bounds...)
text!(ax, 129.4, 38.1, text="E", fontsize=32, font=:bold, color=:grey20)
heatmap!(ax, environmental_layers[1], colormap=[Makie.ColorSchemes.lipari[1]])
heatmap!(ax, total_priority, colormap=:lipari)
poly!(ax, north_korea_polygon, color=:grey85, strokewidth=2, strokecolor=:grey40)
scatter!(ax, bon.coordinates, color=colorant"#fff", strokecolor=:black, strokewidth=2, marker=marker, markersize=markersize)
add_colorbar!(g_top_right[1,1], Makie.ColorSchemes.lipari, title="Total Priority", titlesize=15, height=0.023, width=0.3, halign=0.9, valign=0.22, titlealign=:right)
leg = axislegend(
    ax,
    [MarkerElement(marker = m, color=:white, strokewidth=2, markersize=23) for m in markers],
    ["Occurrence", "Prevalence", "Both"],
    position=:rb,
    "Sampling Type",
    rowgap = 5,
)

f
end

save("plots/korea.png", f)

# ------------------------------------------------------------------------------------
# Make South Korea component of sampling sites figure for box  
# ------------------------------------------------------------------------------------

function get_bbox_of_node(node)
    x,y = SDT.SimpleSDMLayers.__get_grid_coordinate_by_latlon(total_priority, node...)
    Es, Ns = eastings(total_priority), northings(total_priority)

    (Es[x], Es[x+1], Ns[y], Ns[y+1])
end

#node_idx = [3,5,15]
node_idx = [10, 29, 32]
nodes = bon.coordinates[:,node_idx]

for i in node_idx
    @info i,marker[i], get_bbox_of_node(bon.coordinates[:,i])
end

function parse_cities_json(path)
    cities_json = open(path, "r") do f
        return JSON.parse(f)
    end
    
    return [
        Dict(
            :lat => e["lat"],
            :lon => e["lon"],
            :name => e["tags"]["name:en"],
            :population => parse(Int, e["tags"]["population"])
        )     
        for e in cities_json["elements"]]
end 

cities = parse_cities_json(joinpath("data", "korea_cities.json"))
top_num_cities = 10
min_pop = sort([x[:population] for x in cities])[end-top_num_cities+1]

big_cities = filter(x->x[:population] >= min_pop , cities)


#=
    This file is obtained via: 
    
    curl -X POST https://overpass-api.de/api/interpreter -d '
        [out:json][timeout:180];
        area["ISO3166-1"="KR"][admin_level=2]->.sk;
        (
        way(area.sk)["highway"~"^(motorway|trunk|primary)$"];
        way(area.sk)["highway"~"^(motorway_link|trunk_link|primary_link)$"];
        );
        out geom;
        ' > ./data/south_korea_major_highways.json

    and then running 

    `osmtogeojson data/south_korea_major_highways.json > data/south_korea_major_highways.geojson`.
=#

GJ = SDT.GeoJSON
fc = GJ.read(joinpath("data", "south_korea_major_highways.geojson"))


cities_to_plot = ["Seoul", "Busan", "Daegu", "Daejeon", "Gwangju"]

function plot_city_name!(ax, city_name; x_offset = 0.06, y_offset = -0.02)
    idx = findfirst([x[:name] == city_name for x in big_cities])
    (lon, lat) = big_cities[idx][:lon], big_cities[idx][:lat]
    text!(ax, lon+x_offset, lat+y_offset, text=city_name, fontsize=18)
end 

markercols = [colorant"#5390ffff", colorant"#53b37d", colorant"#996bbe"]


@info "\tVisualizing South Korea sites for Box 3 Figure..."
begin 
    f = Figure(size=(700,900))
    ax = Axis(f[1,1], aspect=DataAspect(), xgridvisible=false, ygridvisible=false)
    limits!(ax, bounds...)
    heatmap!(ax, environmental_layers[1], colormap=[:grey90])
    poly!(ax, north_korea_polygon, color=:grey85, strokewidth=2, strokecolor=:grey40)
    lines!(ax, fc.geometry, color=:grey60, linewidth=0.35)
    scatter!(ax, [n for n in nodes], color=markercols, strokecolor=:black, strokewidth=3, marker=marker[node_idx], markersize=30)
    scatter!([(x[:lon], x[:lat]) for x in big_cities], markersize=16, color=:grey30)

    for c in cities_to_plot
        plot_city_name!(ax, c)
    end

    add_colorbar!(g_top_right[1,1], Makie.ColorSchemes.lipari, title="Total Priority", titlesize=15, height=0.023, width=0.35, halign=0.95, valign=0.25, titlealign=:right)
    leg = axislegend(
        ax,
        [MarkerElement(marker = m, color=:white, strokewidth=2, markersize=26) for m in markers],
        ["Occurrence", "Prevalence", "Both"],
        position=:rb,
        "Sampling Type",
        rowgap = 3,
        labelsize=24,
        titlesize=24
    )
    f
end 

save("plots/korea_selected_points.png", f)

[10, 29, 32]

f = Figure()
ax= Axis(f[1,1], aspect=DataAspect())
heatmap!(ax, environmental_layers[1], colormap=[:grey90])
scatter!(ax, bon.coordinates, color=colorant"#fff", strokecolor=:black, strokewidth=1, marker=marker, markersize=8)
for i in 1:size(bon.coordinates, 2)
    text!(ax, bon.coordinates[1,i], bon.coordinates[2,i], text="$i")
end 
f