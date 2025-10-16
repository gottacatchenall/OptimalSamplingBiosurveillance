using KernelDensity
using Distributions
using SpeciesDistributionToolkit
using CairoMakie
using GaussianProcesses
using Optim

const SDT = SpeciesDistributionToolkit
const BONs = BiodiversityObservationNetworks


include("util.jl")



function get_prevalence_data(country_df, host, virus)
    hostvirus_df = filter(x->x.host_species == host && x.pathogen_species_cleaned == virus, country_df)

    nrow(hostvirus_df) > 0 || return nothing
    

    # Get the coordinates of samples
    pts = [(r.longitude, r.latitude) for r in eachrow(hostvirus_df)]
    sample_size = hostvirus_df.number_tested
    # Get prevalence at each sample
    prev = [r.number_positive/r.number_tested for r in eachrow(hostvirus_df)]

    return Dict(
        :host => host,
        :virus => virus,
        :coordinates => pts,
        :sample_size => sample_size,
        :prevalence => prev,
    )
end


function get_prevalence_data(country_df)
    
    hosts = get_hosts(country_df)
    viruses = get_viruses(country_df)

    pairs = []
    for host in hosts
        for virus in viruses
            res = get_prevalence_data(country_df, host, virus)
            if !isnothing(res)
                push!(pairs, res)
            end 
        end
    end 
    return pairs
end


# Get country polygon
countries = getpolygon(PolygonData(NaturalEarth, Countries))
sk_poly = countries["South Korea"]

avg_temp = SDMLayer(RasterData(CHELSA2, AverageTemperature); SDT.boundingbox(sk_poly)...)

sk_df = get_country_df(load_dataframe()
, "South Korea")

hosts = ["Rattus norvegicus", "Apodemus agrarius", "Crocidura lasiura"]
viruses = ["Orthohantavirus hantanense", "Orthohantavirus seoulense"]


hv_pairs = filter(!isnothing, vec([get_prevalence_data(sk_df, h, v) for h in hosts, v in viruses]))




# using GP
pair_idx = 2
lines(sk_poly)
scatter!(hv_pairs[pair_idx][:coordinates], color=hv_pairs[pair_idx][:prevalence])
current_figure()


_euclidean_dist(x,y) = sqrt(sum((x .- y).^2))
function _get_length_scale(coords)
    dist_mat = [_euclidean_dist(coords[i], coords[j]) for i in eachindex(coords), j in eachindex(coords)]
    length_scale = Float64(median(dist_mat))
    return length_scale
end

function fit_gp(coords, prevalence)
    μ = MeanConst(Float64(mean(val)))
    #length_scale = _get_length_scale(coords)
    kern = SE(0., Float64.(log(std(prevalence))))

    coord_mat = Float64.(hcat([[c[1],c[2]] for c in coords]...))

    gp = GP(coord_mat,prevalence,μ,kern)
    optimize!(gp, domean=false)


    full_coords = hcat([[e,n] for e in eastings(avg_temp), n in northings(avg_temp)]...)

    y_predict, σ² = predict_y(gp, full_coords)

    gp_result = Float64.(deepcopy(avg_temp))
    gp_var = Float64.(deepcopy(avg_temp))

    for (i,c) in enumerate(eachcol(full_coords))
        gp_result[c[1],c[2]] = y_predict[i]
        gp_var[c[1],c[2]] = σ²[i]
    end

    mask!(gp_result, sk_poly)
    mask!(gp_var, sk_poly)


    return gp_result, gp_var
end 


hv_pairs = get_prevalence_data(get_country_df(load_dataframe(), "South Korea"))

interpolated_prev = []

for p in hv_pairs
    num_obs = length(p[:prevalence])
    if num_obs > 5 && sum(p[:prevalence]) > 0
        gp_result, gp_var = fit_gp(p[:coordinates], p[:prevalence])
        push!(
            interpolated_prev,
            Dict(
                :host=>p[:host],
                :virus=>p[:virus],
                :prev_predict=>gp_result,
                :prev_var=>gp_var,
                :sampled_coordinates => p[:coordinates]
            )
        )
    end
end



total_prev_uncert = sum([p[:prev_var] for p in interpolated_prev])
heatmap(total_prev_uncert)

idx = 4
heatmap(interpolated_prev[idx][:prev_var])
scatter!(interpolated_prev[idx][:sampled_coordinates])
current_figure()



# combine kriging w/ sdm

korea_sdms = read_sdms("SDM_Artifacts_SouthKorea")

uncs = [quantize(v[:uncertainty]) for (k,v) in korea_sdms]
pred = [quantize(v[:prediction]) for (k,v) in korea_sdms]


species_weights = rand(length(uncs))
species_weights ./= sum(species_weights)

sdm_unc = sum([species_weights[i]*uncs[i] for i in eachindex(uncs)])
sdm_pred = sum([species_weights[i]*pred[i] for i in eachindex(pred)])


total_prev_uncert.x = sdm_unc.x
total_prev_uncert.y = sdm_unc.y


hist(quantize(sdm_unc) + quantize(sdm_pred) + quantize(total_prev_uncert))



w_sdm = 0.5 

priority = w_sdm*0.5*(quantize(sdm_unc) + quantize(sdm_pred) )  + (1-w_sdm) * quantize(total_prev_uncert)

extrema(priority)

hist(priority)
extrema(priority)

titled_priority = BONs.tilt(priority, 5.)
hist(titled_priority)

bon = BONs.sample(BalancedAcceptance(), priority, inclusion=titled_priority)

heatmap(titled_priority)
scatter!([n for n in bon], color=:red)
current_figure()

