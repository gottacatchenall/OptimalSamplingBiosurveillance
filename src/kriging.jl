using Distributions
using SpeciesDistributionToolkit
using CairoMakie
using GaussianProcesses
using Optim

const SDT = SpeciesDistributionToolkit

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





