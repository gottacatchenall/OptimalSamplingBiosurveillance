using CSV, DataFrames, HTTP, JSON
using CairoMakie
using SpeciesDistributionToolkit
using Statistics

using TernaryDiagrams

using BiodiversityObservationNetworks

const SDT = SpeciesDistributionToolkit

include("util.jl")
include("gbif.jl")

function split_occurrences_into_species(occs, min_occs=50)
    function _transform_to_genus_species(str)
        x = split(str, " ")[1:2]
        return x[1]*" "*x[2]
    end
    unique_species = unique(map(
        x->_transform_to_genus_species(x.what),
        occs)
    )

    dict = Dict()
    for sp in unique_species
        these_occs = Occurrences(filter(x->_transform_to_genus_species(x.what) == sp, SDT.elements(occs)))
        @info sp, ": ", length(these_occs)
        if length(these_occs) > min_occs
            dict[sp] = these_occs
        end
    end
    return dict
end

function get_pseudoabsences(
    presence_layer, 
    buffer_distance,
    class_balance
)
    background = pseudoabsencemask(DistanceToEvent, presence_layer)
    
    buffer = pseudoabsencemask(WithinRadius, presence_layer; distance = buffer_distance)

    ndmask = nodata(buffer, true) 
    background.indices .= ndmask.indices

    pseudoabs_layer = backgroundpoints(background, Int(round(class_balance*sum(presence_layer))))
    return pseudoabs_layer
end

function validation_metrics(confusion_matrices)
    return Dict(
        :mcc => mean(mcc.(confusion_matrices)),
        :tss => mean(trueskill.(confusion_matrices)),
        :κ => mean(κ.(confusion_matrices))
    )
end

function fit_sdm(
    occurrences,
    layers;
    pseudoabsence_buffer_distance = 5.,
    class_balance = 1., # high val -> more negatives
    num_trees = 10,
    num_folds = 3,
)
    @info "Starting SDM for $(occurrences[1].what)"

    presence_layer = mask(layers[begin], occurrences)
    absence_layer = get_pseudoabsences(presence_layer, pseudoabsence_buffer_distance, class_balance)



    #dt = SDM(RawData, DecisionTree, layers, presence_layer, absence_layer)
    dt = SDM(SDT.ZScore, SDT.Logistic, layers, presence_layer, absence_layer)

    #random_forest = Bagging(dt, num_trees)
    #bagfeatures!(random_forest)

    # var selection
    variables!(dt, ForwardSelection, kfold(dt, k=num_folds))

    train!(dt)

    cv_results = crossvalidate(dt, kfold(dt, k=num_folds))

    unc_prs = [partialresponse(dt, L, rand(variables(dt)), threshold=false, inflated=true) for _ in 1:100]


    unc = mosaic(var, unc_prs)

    return Dict(
        :model => dt,
        :presences => presence_layer,
        :absences => absence_layer,
        :range => SDT.predict(dt, layers, threshold=true),
        :prediction => SDT.predict(dt, layers, threshold=false),
        :uncertainty => unc,
        :confusion => cv_results,
        :metrics => validation_metrics(cv_results.validation)
    )
end

function make_sdms(
    occs,
    layers;
    min_occs = 50,
    kwargs...
)
    occs_per_species = split_occurrences_into_species(occs, min_occs)
    
    return Dict([k=>fit_sdm(v, layers; kwargs...) for (k,v) in occs_per_species])
end



"""
load_env()
COUNTRY = "South Korea"
full_df = load_dataframe()
country_df = get_country_df(full_df, COUNTRY)

hosts = get_hosts(country_df)
viruses = get_viruses(country_df)


# TODO figure out the GBIF.jl way here
#hacky_get_gbif_data(hosts, country_df)
"""

# TODO: have the doi yanked from the GBIF download center
doi = "10.15468/dl.3g2ppt"
occ = GBIF.download(doi)

countr = getpolygon(PolygonData(NaturalEarth, Countries))
sk = countr["South Korea"]


L = [Float32.(SDMLayer(RasterData(CHELSA2, BioClim); layer=i, SDT.boundingbox(sk)...)) for i in 1:19]
mask!(L, sk)

@time res = make_sdms(
    occ, 
    L;
    pseudoabsence_buffer_distance = 20.,
    class_balance = 2., # high val -> more negatives
    num_trees = 30,
    num_folds = 10
)


write_sdm_artifacts("SDM_Artifacts_SouthKorea", res)


"""
# --------------------------------------------------------
# India SDMs


doi = "10.15468/dl.b3z33h"
occ = GBIF.download(doi)



countr = getpolygon(PolygonData(NaturalEarth, Countries))
india_poly = countr["India"]

L = [Float32.(SDMLayer(RasterData(WorldClim2, BioClim); resolution=2.5, layer=i, SDT.boundingbox(india_poly)...)) for i in 1:19]

mask!(L, india_poly)



@time res = make_sdms(
    occ, 
    L;
    pseudoabsence_buffer_distance = 20.,
    class_balance = 2., # high val -> more negatives
    num_trees = 30,
    num_folds = 10,
    min_occs = 25
)



w = rand(length(res))
w ./= sum(w)

total_unc = sum([w[i]*quantize(r[:uncertainty]) for (i,r) in enumerate(values(res))])

total_pred = sum([w[i]*quantize(r[:prediction]) for (i,r) in enumerate(values(res))])

scatter(total_pred, total_unc, color=(:blue, 0.01))


bon = BiodiversityObservationNetworks.sample(BalancedAcceptance(), total_pred, inclusion=BiodiversityObservationNetworks.tilt(total_pred+total_unc, 3))


hist(total_unc + total_pred)

heatmap(total_unc + total_pred)
scatter!([n for n in bon.nodes], color=:red)
current_figure()

write_sdm_artifacts("SDM_Artifacts_India", res)


occs_per_species = split_occurrences_into_species(occ)


occs = occs_per_species["Rattus rattus"]

res = fit_sdm(
    occs, L;
    pseudoabsence_buffer_distance = 20.,
    class_balance = 2., # high val -> more negatives
    num_trees = 30,
    num_folds = 10,
)

model = res[:model]
train!(model)

heatmap(SDT.predict(model, L, threshold=false))
scatter!(res[:presences], color=:red)
current_figure()


prs = [partialresponse(model, L, rand(variables(model)), threshold=false, inflated=true) for _ in 1:100]


heatmap(mosaic(var, prs))
scatter!(res[:presences], color=:red)
current_figure()

scatter(quantize(res[:prediction]), quantize(mosaic(var, prs)), color=(:blue, 0.03))



partialresponse(model, L, variables(model)[4], threshold=false) |> heatmap


heatmap(res[:prediction])







res = read_sdms("SDM_Artifacts_India")
sp = "Funambulus tristriatus"

heatmap(res[sp][:prediction])
scatter!(Bool.(res[sp][:presences]), color=:red)
current_figure()




countr = getpolygon(PolygonData(NaturalEarth, Countries))
india_poly = countr["India"]

lines(india_poly)

"""
