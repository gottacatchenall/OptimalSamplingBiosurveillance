using CSV, DataFrames, JSON

function write_sdm_artifacts(artifact_dir, results)
    mkpath(artifact_dir)

    for (sp_name,v) in results
        output_dir = joinpath(artifact_dir, sp_name)
        mkpath(output_dir)

        SDT.SimpleSDMLayers.save(
            joinpath(output_dir, "presences.tif"),
            Float32.(v[:presences])
        )      
        SDT.SimpleSDMLayers.save(
            joinpath(output_dir, "absences.tif"),
            Float32.(v[:absences])
        )      
        SDT.SimpleSDMLayers.save(
            joinpath(output_dir, "uncertainty.tif"),
            v[:uncertainty]
        )        
        SDT.SimpleSDMLayers.save(
            joinpath(output_dir, "prediction.tif"),
            v[:prediction]
        )        

        open(joinpath(output_dir, "metrics.json"), "w") do f
            JSON.print(f, v[:metrics])
        end

    end

end

# Read the SDMs
function read_sdms(dir)
    species_names = readdir(dir)
    Dict([
        s=> Dict(
                :uncertainty => SDMLayer(joinpath(joinpath(dir, s), "uncertainty.tif")),
                :prediction => SDMLayer(joinpath(joinpath(dir, s), "prediction.tif")),
                :metrics => open(joinpath(joinpath(dir, s), "metrics.json"), "r") do f
                    return JSON.parse(f)
                end,
                :presences => SDMLayer(joinpath(joinpath(dir, s), "presences.tif")),
                :absences => SDMLayer(joinpath(joinpath(dir, s), "absences.tif")),
            ) 
            for s in species_names
        ]
    )
end

function get_prevalence_data()
    df = CSV.read(joinpath("data", "prevalence.csv"), DataFrame)
    hosts, viruses = unique(df.host), unique(df.virus)
    result = []
    for h in hosts, v in viruses
        this_pair_df = filter(x->x.host == h && x.virus == v, df)
        if nrow(this_pair_df) > 0 
            push!(
                result, 
                Dict(
                    :host => h,
                    :virus => v,
                    :coordinates => [(this_pair_df.longitude[i], this_pair_df.latitude[i]) for i in eachindex(this_pair_df.prevalence)],
                    :prevalence => this_pair_df.prevalence,
                )
            )
        end 
    end
    return result
end
