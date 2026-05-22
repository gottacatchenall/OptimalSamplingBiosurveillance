# Optimal Sampling for Biosurveillance

This repository contains code associated with the manuscript _A protocol for biodiversity-informed wildlife disease surveillance_ 

# Running the Code

This analysis was written in Julia v1.12. Instructions for installing Julia can be found [here](https://julialang.org/downloads/).  

The required packages can be installed by running Julia from the root of the repository, and running

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Before running all analysis, make sure to follow the instructions in the below section
(Obtaining South Korea highway data), to ensure the visualization in Box 3 can ve built
properly.

All analysis can then be run with `julia --project=. main.jl`

# Obtaining South Korea highway data

For the figure in Box 3, the South Korea highway network data is too large to store in the repo. It can be obtained via running to query OpenStreetMap

```
curl -X POST https://overpass.kumi.systems/api/interpreter \
        --data-urlencode 'data=
          [out:json][timeout:180];
          area["ISO3166-1"="KR"][admin_level=2]->.sk;
          (
          way(area.sk)["highway"~"^(motorway|trunk|primary)$"];
          way(area.sk)["highway"~"^(motorway_link|trunk_link|primary_link)$"];
          );
          out geom;
        ' > ./data/south_korea_major_highways.json
```

and then running 

```
osmtogeojson ./data/south_korea_major_highways.json > ./data/south_korea_major_highways.geojson
```

to convert to GeoJSON.

`osmtogeojson` will need to be installed with `npm install -g osmtogeojson`. If you don't have `npm`, instructions for installation can be found [here](https://docs.npmjs.com/downloading-and-installing-node-js-and-npm).

# File Structure

## `./src`

- `sdms.jl`: Contains functions for fitting and validating species distribution models
- `kriging.jl`: Contains functions for Gaussian Process regression (aka kriging) for estimating prevalence across space for the South Korea case study
- `util.jl`: Utilities for reading and writing SDMs, and querying the data files
- `plotting.jl`: Utilities for plotting 


## `./scripts`

- `india.jl`: Code to make the case study for India
- `southkorea.jl`: Code to make the case study for South Korea
- `sensitivity_analysis.jl`: Code to run sensitivity analysis
- `sensitivity_visualization.jl`: Code to make the visualization for sensitivity analysis


## `./data`:
- `lassa_india.csv`: The survey data from Rodrigues et al. (1978) testing for Lassa virus in India
- `prevalence.csv`: Prevalence data for the South Korea Case Study
- `korea_cities.json`: City locations for the Figure in Box 3


## `./osm`:
- `osm_viz.py`: A python script to obtain OpenStreetMap data for the sampling locations, for use in making the Box 3 figure.
