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




function load_dataframe(
    host_path = joinpath("data", "host_clean.csv"),
    pathogen_path = joinpath("data", "pathogen_clean.csv"),
    columns_to_select = [
        "host_record_id",
        "study_id",
        "host_species",
        "country",
        "latitude",
        "longitude",
        "gbif_id",
        "coordinate_resolution_processed",
        "iso3c",
        "start_date",
        "end_date",
        "pathogen_species_cleaned",
        "number_tested",
        "number_positive"
    ]
)

    columns_to_convert_to_float = [
        "latitude",
        "longitude",
        "number_tested",
        "number_positive"
    ]

    host_df = CSV.read(host_path, DataFrame)
    pathogen_df = CSV.read(pathogen_path, DataFrame)
    host_df, pathogen_df

    df = leftjoin(
        host_df, 
        pathogen_df, 
        on=:host_record_id=>:host_record_id, 
        makeunique=true
    )

    for c in columns_to_select
        filter!(row->!ismissing(row[c]) && row[c] != "NA", df)
    end 
    #filter!(row->!ismissing(row.pathogen_species_cleaned) && row.pathogen_species_cleaned != "NA" && row.host_species != "NA" && row.latitude != "NA" && row.longitude != "NA", df)
    
    iso_df = CSV.read(joinpath("data", "iso_codes.csv"), DataFrame)

    df = select(df, columns_to_select)

    function _match_iso3d(iso3)
        iso3 == "NA" && return missing
        idx = findfirst(r->r["Alpha-3 code"] == iso3, eachrow(iso_df))
        return iso_df[idx,"Alpha-2 code"]
    end
    df.iso2d = _match_iso3d.(df.iso3c)

    for c in columns_to_convert_to_float
        df[!,c] .= parse.(Float32, df[!,c])
    end

    return df
end

function get_country_df(full_df, country_name)
    country_name ∈ full_df.country || throw(ArgumentError("Country $country_name is not in the DataFrame"))

    country_df = filter(row->row.country == country_name, full_df)
    # Drop all locations without tests
    filter!(x->x.number_tested > 0, country_df)

    # Drop all samples where the sk_df is coarse (country-level)
    filter!(x->x.coordinate_resolution_processed != "country", country_df)
    return country_df
end

function get_hosts(df)
    String.(unique(df.host_species))
end

function get_hosts(df, virus)
    this_virus_df = filter(row->Base.contains(row.pathogen_species_cleaned, virus), df)

    unique(this_virus_df.host_species)
end



function get_viruses(df)
    unique_lines = unique(df.pathogen_species_cleaned)
    String.(unique(vcat(split.(unique_lines, ", ")...)))
end

function get_viruses(df, host)
    this_host_df = filter(row->row.host_species == host, df)
    return unique(this_host_df.pathogen_species_cleaned)
end



df = load_dataframe()

filter(c->c.country=="Senegal", df).pathogen_species_cleaned

df = get_country_df(df, "Senegal")

df

df.number_positive

df.pathogen_species_cleaned