using CSV, DataFrames, HTTP, JSON
using CairoMakie
using SpeciesDistributionToolkit
using Statistics
using ColorSchemes
using TernaryDiagrams
using CairoMakie.Colors

using BiodiversityObservationNetworks

const SDT = SpeciesDistributionToolkit

include("util.jl")


korea_sdms = read_sdms("SDM_Artifacts_SouthKorea")

richness = sum([v[:prediction] for (k,v) in korea_sdms])
uncert = sum([v[:uncertainty] for (k,v) in korea_sdms])


function get_quantiles(layer; num_bins=5)
    quantile(layer, [x for x in LinRange(0,1, num_bins+1)][2:end-1])         
end

function quantize_layer(layer, quantiles)
    newlayer = similar(layer)
    for i in eachindex(layer.grid)
        tmp = findfirst(x-> layer.grid[i] < x, quantiles)
        newlayer.grid[i] = !isnothing(tmp) ? tmp : length(quantiles)+1
    end
    newlayer 
end

num_bins = 8

quant_richness, quant_uncert = get_quantiles(richness; num_bins=num_bins), get_quantiles(uncert; num_bins=num_bins)


richness_quantized, uncert_quantized = quantize_layer(richness, quant_richness), quantize_layer(uncert, quant_uncert)

encoding = ([num_bins*j+i for i in 1:num_bins, j in 0:(num_bins-1)])

bivar_layer = similar(richness)
for I in eachindex(uncert)
    bivar_layer.grid[I] = encoding[Int32(richness_quantized.grid[I]), Int32(uncert_quantized.grid[I]),]
end
bivar_layer


function setup_axis(ax, layer)
    EXTENT = SDT.boundingbox(layer)

    long = EXTENT[:left], EXTENT[:right]
    lat = EXTENT[:bottom], EXTENT[:top]

    blue_red_4 = vec(reshape([color(x) for x in ["#d3d3d3", "#cda2a2", "#c66c6e", "#be2427", "#adc3cb", "#a8959c", "#a2646a", "#9c2226", "#87b2c3", "#838995", "#7f5b66", "#791f24", "#5fa1ba", "#5c7c8e", "#595361", "#551c22"]], (4,4)))
    blue_red_8 = [color(x) for x in ["#d3d3d3", "#cfbcbd", "#c9a6a6", "#c59091", "#c07879", "#bb6061", "#b44446", "#ad2123", "#c2cacd", "#beb4b8", "#b99fa1", "#b48a8c", "#b07375", "#ab5b5e", "#a54144", "#9f1f22", "#b0c1c7", "#acacb3", "#a7989d", "#a38389", "#a06e72", "#9b575c", "#963e42", "#901e21", "#9eb8c1", "#9ba4ad", "#979198", "#937d84", "#90696f", "#8c5359", "#873b40", "#821d20", "#8cafbb", "#899ca8", "#858a94", "#827780", "#7f636b", "#7c4f56", "#78383e", "#731b20", "#7ba5b5", "#7893a2", "#75828f", "#72707c", "#705e68", "#6d4b54", "#69353c", "#651a1e", "#699cb0", "#668b9d", "#647b8b", "#616a78", "#5f5965", "#5c4651", "#59323a", "#56181e", "#5692a9", "#558398", "#527485", "#506474", "#4f5361", "#4c424e", "#4a2f38", "#47171c"]]
    blue_red_5 = vec(reshape([color(x) for x in ["#d3d3d3", "#cbacad", "#c28485", "#b9595b", "#ad2123", "#b4c3c9", "#ad9fa5", "#a57a7f", "#9e5257", "#941e22", "#95b3be", "#8f929c", "#897078", "#834c52", "#7a1c20", "#76a3b4", "#728594", "#6d6671", "#67454e", "#61191e", "#5692a9", "#53778b", "#4f5c6a", "#4c3e49", "#47171c"]], (5,5)))


    hm = heatmap!(ax, layer, colormap=blue_red_8)
end 



blue_red_8 = [color(x) for x in ["#d3d3d3", "#cfbcbd", "#c9a6a6", "#c59091", "#c07879", "#bb6061", "#b44446", "#ad2123", "#c2cacd", "#beb4b8", "#b99fa1", "#b48a8c", "#b07375", "#ab5b5e", "#a54144", "#9f1f22", "#b0c1c7", "#acacb3", "#a7989d", "#a38389", "#a06e72", "#9b575c", "#963e42", "#901e21", "#9eb8c1", "#9ba4ad", "#979198", "#937d84", "#90696f", "#8c5359", "#873b40", "#821d20", "#8cafbb", "#899ca8", "#858a94", "#827780", "#7f636b", "#7c4f56", "#78383e", "#731b20", "#7ba5b5", "#7893a2", "#75828f", "#72707c", "#705e68", "#6d4b54", "#69353c", "#651a1e", "#699cb0", "#668b9d", "#647b8b", "#616a78", "#5f5965", "#5c4651", "#59323a", "#56181e", "#5692a9", "#558398", "#527485", "#506474", "#4f5361", "#4c424e", "#4a2f38", "#47171c"]]

f = Figure(resolution=(2000,2000))
ax = Axis(f[1,1],
    aspect=DataAspect(),
    titlealign=:left,
    yticklabelsvisible=false,
    xticklabelsvisible=false,
    xgridvisible=false,
    ygridvisible=false
)
setup_axis(ax, bivar_layer)
#scatter!(ax, [n for n in bon], color=:orange, markersize=35)

ax_inset = Axis(f[1, 1],
width=Relative(0.2),
height=Relative(0.2),
aspect=1,
halign=0.9,
valign=0.3,
)


heatmap!(ax_inset, (reshape(blue_red_8, (8,8))))
hidedecorations!(ax)

f   

save("bivar.png", f)