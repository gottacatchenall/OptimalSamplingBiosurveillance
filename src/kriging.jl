using Distributions
using SpeciesDistributionToolkit
using CairoMakie
using GaussianProcesses
using Optim

const SDT = SpeciesDistributionToolkit

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
        :study_id=>hostvirus_df.study_id
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





_euclidean_dist(x,y) = sqrt(sum((x .- y).^2))
function _get_length_scale(coords)
    dist_mat = [_euclidean_dist(coords[i], coords[j]) for i in eachindex(coords), j in eachindex(coords)]
    length_scale = Float64(median(dist_mat))
    return length_scale
end

function fit_gp(coords, prevalence, template, polygon)
    μ = MeanConst(Float64(mean(prevalence)))
    kern = SE(0., Float64.(log(std(prevalence))))

    coord_mat = Float64.(hcat([[c[1],c[2]] for c in coords]...))

    gp = GP(coord_mat,prevalence,μ,kern)
    optimize!(gp, domean=false)


    full_coords = hcat([[e,n] for e in eastings(template), n in northings(template)]...)

    y_predict, σ² = predict_y(gp, full_coords)

    gp_result = Float64.(similar(template))
    gp_var = Float64.(similar(template))

    for (i,c) in enumerate(eachcol(full_coords))
        gp_result[c[1],c[2]] = y_predict[i]
        gp_var[c[1],c[2]] = σ²[i]
    end

    mask!(gp_result, polygon)
    mask!(gp_var, polygon)


    return gp_result, gp_var
end 





