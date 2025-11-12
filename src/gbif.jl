function load_env()
    dotenv = open(".env", "r") do file
        read(file, String)
    end
    envvars = split(dotenv, "\n")
    
    for ev in envvars
        k, v = split(ev, "=")
        ENV[k] = v
    end
end


function build_taxon_query(taxa)
    ts = GBIF.taxon.(taxa)
    taxon_query = []
    for t in ts
        levels = [:kingdom, :phylum, :class, :order, :family, :genus, :species]
        level = levels[findlast(l -> getfield(t, l) !== missing, levels)]
        push!(taxon_query, String(level) * "Key" => getfield(t, level).second)
    end
    return taxon_query
end

function get_gbif_data(country_df, host_species, notification=true)
    country_code = country_df.iso2d[begin]
    tx_query = build_taxon_query(host_species)
    query =
        [
            tx_query...,
            "hasCoordinate" => true,
            "country" => country_code,
            "occurrenceStatus" => "PRESENT",
        ]

end 


function hacky_get_gbif_data(taxa, country_code)
    #tx_query = build_taxon_query(host_species)
    ts = GBIF.taxon.(taxa)
    species_ids = [t.species[2] for t in ts]

    query_template = JSON.parse(String(read(joinpath("data", "gbif_query_template.json"))))

    query_template["creator"] = ENV["GBIF_USERNAME"]
    
    push!(
        query_template["predicate"]["predicates"], Dict([
            "type" => "in",
            "key" => "COUNTRY",
            "values" => [country_code]
        ])
    )

    push!(
        query_template["predicate"]["predicates"], Dict([
            "type" => "in",
            "key" => "TAXON_KEY",
            "values" => species_ids
        ])
    )
    usn, pwd = ENV["GBIF_USERNAME"], ENV["GBIF_PASSWORD"]

    query_path = joinpath(tempdir(), "query.json")


    open(query_path, "w") do io
        JSON.print(io, query_template)
    end 

    cmd = `curl --include --user $(usn):$(pwd) --header "Content-Type: application/json" --data @$(query_path) https://api.gbif.org/v1/occurrence/download/request`
    run(cmd)
end
