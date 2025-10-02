using SpeciesDistributionToolkit
using CSV, DataFrames
using CairoMakie
using BiodiversityObservationNetworks
using Statistics

const BONs = BiodiversityObservationNetworks
const SDT = SpeciesDistributionToolkit


doi = "10.15468/dl.q3wd34"
occ = GBIF.download(doi)

species = ["Mus musculus", "Myodes regulus", "Crocidura lasiura", "Apodemus agrarius"]

# Remove records without coordinates
filter!(x -> !ismissing(place(x)), elements(occ))

# Remove the point in the pacific
filter!(x -> longitudes(x) > 100, elements(occ))


occurrence_records =
    [Occurrences(filter(x -> startswith(s)(entity(x)), elements(occ))) for s in species]

# Get a polygon for South Korea
countr = getpolygon(PolygonData(NaturalEarth, Countries))
sk = countr["South Korea"]

# Visualization of the first species occurrences
poly(sk)
scatter!(occurrence_records[begin], color = :red, markersize = 3)
current_figure()


# Download CHELSA2 BioClim Layers
L = [
    Float32.(SDMLayer(RasterData(CHELSA2, BioClim)); layer = i, SDT.boundingbox(sk)...)
    for i = 1:19
]

# Mask layers to South Korea's Land
mask!(L, sk)

# Convert occurrence records to layers
occ_layer = [mask(L[1], o) for o in occ_dict]

# Fits an SDM for the species at index `sp_idx` in the `species` array 
function fit_sdm(sp_idx)
    occ_i = occ_layer[sp_idx]

    pseudoabs_surface = pseudoabsencemask(DistanceToEvent, occ_i)
    bgpoints = backgroundpoints(pseudoabs_surface, sum(occ_i))

    sdm = SDM(RawData, NaiveBayes, L, occ_i, bgpoints)
    bag_sdm = Bagging(sdm, 30)
    train!(bag_sdm)
    return bag_sdm
end

# Fit SDMs (2 is omitted due to lack of presences)
sdms = fit_sdm.([1, 3, 4])


# lightweight variable importance metric
vimp = variableimportance(sdms[1], kfold(sdms[1], k = 2))


# Weights uncertainty with more focus toward high values (i.e. species is likely present)
# (This is no longer used)
weight_uncertainty(x) = (minimum(x) + 4 * Statistics.median(x) + maximum(x)) / 6

# Predict SDMs
pred = [predict(sdm, L) for sdm in sdms]

# Compute uncertainty as bootstrap variance (measured as the range between the 25% and 75% quantiles, iqr)
unc = [predict(sdm, L, consensus = iqr) for sdm in sdms]
#unc = [predict(sdm, L, consensus=weight_uncertainty) for sdm in sdms]


# Randomly weight importance of each species
r = rand(3)
r = r ./ sum(r)
weighted_pred = [r[i] * pred[i] for i in eachindex(pred)]
weighted_unc = [r[i] * unc[i] for i in eachindex(unc)]


# Compute the total prediction score
total_pred = sum(weighted_pred)

# Compute the total uncertainty
total_unc = mosaic(var, weighted_unc)


# Split into quadrants 
p, q = 0.5, 0.5
pred_cutoff = quantile(total_pred, [p])[begin]
unc_cutoff = quantile(total_unc, [q])[begin]

A = total_pred .> pred_cutoff
B = total_unc .> unc_cutoff
bon = sample(BalancedAcceptance(), total_unc, inclusion = total_unc)

heatmap(nodata(A & B, false), colormap = ["#b48ead"], label = "Active Surv.")
heatmap!(nodata(.!A & .!B, false), colormap = ["#5e81ac"], label = "Confirmation")
heatmap!(nodata(.!A & B, false), colormap = ["#a3be8c"], label = "Discovery")
heatmap!(nodata(A & .!B, false), colormap = ["#d08770"], label = "Passive Surv.")
axislegend()
scatter!([n for n in bon.nodes], color = :black)

current_figure()
