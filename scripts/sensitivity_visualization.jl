@info "Making Sensitivity Visualization..."

using ColorBlendModes
using Random

include(joinpath("..", "src", "sdms.jl"))
include(joinpath("..", "src", "plotting.jl"))
include(joinpath("..", "src", "util.jl"))

using BiodiversityObservationNetworks
const AG = SDT.SimpleSDMPolygons.AG

sensitivity_df = CSV.read(joinpath("artifacts", "sensitivity.csv"), DataFrame)

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


# Remove contested areas from India polygon
kashmir_polygon = india_polygon[2]
arunuachal_pradesh_polygon = india_polygon[5]

india_polygon = india_polygon - arunuachal_pradesh_polygon



# Remove contested are
kashmir_polygon = india_polygon[2]
arunuachal_pradesh_polygon = india_polygon[5]
india_polygon = india_polygon - kashmir_polygon - arunuachal_pradesh_polygon

# Get polygons with dashes for contested regions
kashmir_dashes = add_dashes(kashmir_polygon)
ap_dashes = add_dashes(arunuachal_pradesh_polygon)


sdms = read_sdms(joinpath("artifacts", "India"))

# Species to use
species = ["Rattus rattus", "Suncus murinus", "Hystrix indica"]

# Get uncertainties for each species
uncs = [sdms[s][:uncertainty] for s in species]
richness = [sdms[s][:prediction] for s in species]



const r1 = [0, 0]
const r2 = [1, 0]
const r3 = [0.5, sqrt(3) / 2]
const R = [
    1 1 1
    r1 r2 r3
]

function compute_ternary_pt((x,y,z))
    cartidx = R * [x, y, z]
    return (cartidx[2], cartidx[3])
end 

# Adapted from TernaryDiagrams.jl
function draw_grid!(ax::Axis; grid_line_width=0.5, grid_line_color=:grey50)
    # draw grid
    fracs = 0.0:0.1:1.0

    for f1 in fracs
        f2 = 1 - f1
        vec1 = [f1, f2, 0]
        vec2 = [f1, 0, f2]

        x1 = Point2f((R*vec1)[2:3]...)
        x2 = Point2f((R*vec2)[2:3]...)

        lines!(ax, [x1, x2], linewidth = grid_line_width, color = grid_line_color)

        # labelx
        f1 in 0:0.2:1.0 && continue
        text!(
            ax,
            x1,
            text = "  " * string(round(f1, digits = 2)),
            rotation = -π / 3,
            align = (:left, :center),
        )
    end

    for f1 in fracs
        f2 = 1 - f1
        vec1 = [0, f2, f1]
        vec2 = [f1, f2, 0]

        x1 = Point2f((R*vec1)[2:3]...)
        x2 = Point2f((R*vec2)[2:3]...)

        lines!(ax, [x1, x2], linewidth = grid_line_width, color = grid_line_color)

        #labely
        f1 in 0:0.2:1.0 && continue
        text!(
            ax,
            x1,
            text = "  " * string(round(f2, digits = 2)),
            rotation = π / 3,
            align = (:left, :center),
        )
    end

    for f1 in fracs
        f2 = 1 - f1
        vec1 = [f2, 0, f1]
        vec2 = [0, f2, f1]

        x1 = Point2f((R*vec1)[2:3]...)
        x2 = Point2f((R*vec2)[2:3]...)

        lines!(ax, [x1, x2], linewidth = grid_line_width, color = grid_line_color)

        # labelz
        f1 in 0:0.2:1.0 && continue
        text!(
            ax,
            x1,
            text = string(round(f1, digits = 2)) * "  ",
            align = (:right, :center),
        )
    end
end



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

ax_args = (
    xgridvisible=false,
    ygridvisible=false,
    aspect=DataAspect()
)


m, M = extrema(sensitivity_df.mae)
relative_MAE = (sensitivity_df.mae .- m) ./ (M - m)

_3d_pts = [(x,y,z) for (_, x,y,z) in eachrow(sensitivity_df)]
_2d_pts = compute_ternary_pt.(_3d_pts)



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
        [0.1, 0.1, 0.8]
    ]


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

    ax_ternary = Axis(g[1,1], aspect=1)
    ylims!(ax_ternary, -0.1, 0.9)
    hidedecorations!(ax_ternary)
    hidespines!(ax_ternary)
    scatter!(ax_ternary, _2d_pts, color=relative_MAE, marker=:hexagon, markersize=25, colormap=:magma)
    draw_grid!(ax_ternary, grid_line_color=:white, grid_line_width=0.25)

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
    

    text!(ax_ternary, 0.06, 0.32, text="Rattus rattus weight", rotation=π/3, fontsize=20)
    text!(ax_ternary, 0.32, -0.15, text="Suncus murinus weight", fontsize=20)
    text!(ax_ternary, 0.76, 0.65, text="Hystrix indica weight", rotation=5π/3, fontsize=20)
    colsize!(g, 1, Relative(0.7))

    fig

end 

mkpath("plots")
save(joinpath("plots", "sensitivity.png"), fig)
