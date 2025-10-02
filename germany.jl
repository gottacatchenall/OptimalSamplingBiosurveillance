using SpeciesDistributionToolkit
using CSV, DataFrames
using CairoMakie
using BiodiversityObservationNetworks
const BONs = BiodiversityObservationNetworks
const SDT = SpeciesDistributionToolkit

using Statistics

# 
host_df = CSV.read(joinpath("data", "host_clean.csv"), DataFrame)
path_df = CSV.read(joinpath("data", "pathogen_clean.csv"), DataFrame)

full_df =
    leftjoin(host_df, path_df, on = :host_record_id=>:host_record_id, makeunique = true)

# Select only Germany and cast useful columns to Floats
germany_df = filter(x->x.country == "Germany", full_df)
germany_df.number_tested = parse.(Float32, germany_df.number_tested)
germany_df.number_positive = parse.(Float32, germany_df.number_positive)
germany_df.latitude = parse.(Float32, germany_df.latitude)
germany_df.longitude = parse.(Float32, germany_df.longitude)

# Drop all locations without tests
filter!(x->x.number_tested > 0, germany_df)
# Drop all samples where the coordinate is coarse (country-level)
filter!(x->x.coordinate_resolution_processed != "country", germany_df)

# Get a polygon for Germany
countr = getpolygon(PolygonData(NaturalEarth, Countries))
ger = countr["Germany"]

# Get the coordinates of samples
pts = [(r.longitude, r.latitude) for r in eachrow(germany_df)]

# Get prevalence at each sample
prev = [r.number_positive/r.number_tested for r in eachrow(germany_df)]


# Pick a given host
host_idx = findall(x->x.host_species == "Myodes glareolus", eachrow(germany_df))

poly(ger)
scatter!(
    pts,
    color = prev,
    colorrange = extrema(prev),
    markersize = 3log.(germany_df.number_tested .+ 1),
)
current_figure()


# Make a shitty SDM

tx = GBIF.taxon("Myodes glareolus")
query = [
    "hasCoordinate" => true,
    "country" => "DE",
    "limit" => 300,
    "occurrenceStatus" => "PRESENT",
]
occ = GBIF.occurrences(tx, query...)


L = [
    Float32.(SDMLayer(RasterData(CHELSA2, BioClim); layer = i, SDT.boundingbox(ger)...)) for i = 1:19
]
mask!(L, ger)

occ_layer = mask(L[1], occ)
pabs = pseudoabsencemask(DistanceToEvent, occ_layer)
bgpoints = backgroundpoints(pabs, sum(occ_layer))
sdm = SDM(RawData, NaiveBayes, L, occ_layer, bgpoints)
bag_sdm = Bagging(sdm, 5)
train!(bag_sdm)

pred = predict(bag_sdm, L, threshold = false)
unc = predict(bag_sdm, L; consensus = iqr, threshold = false)

heatmap(unc)
scatter!(
    pts[host_idx],
    color = prev[host_idx] .== 0,
    colormap = [:grey85, :black],
    markersize = 15,
)
current_figure()


these_pts = pts[host_idx]

# Compute a buffer distance around each point
dist_mat = [
    SDT.Fauxcurrences._distancefunction(these_pts[i], these_pts[j]) for
    i in eachindex(these_pts), j in eachindex(these_pts)
]
median(dist_mat)
