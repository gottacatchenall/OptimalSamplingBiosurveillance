#import "@preview/colorful-boxes:1.4.3": *
#set text(lang: "en", size:10pt)


#colorbox(
  title: text(font: "Helvetica Neue", weight: 500)[Box 1: Guidelines for the interpretation of Species Distribution Models], 
  color: (
    fill: rgb("#f4fcfa"),
    stroke: rgb("#418b7f"),
    title: rgb("#101010")
  ),
  radius: 2pt,
)[  
  
Species distribution models (SDMs) are widely used tools in ecology that relate species occurrences with environmental variables to identify how species are geographically constrained by different ecological factors @Elith2011StatisticalExplanation. SDMs can then be used to estimate the realized ecological niche of a species --- the environmental space (rather than geographic space) that is suitable for the species to occupy. Because occurrence data is often presence-only, meaning there are few/no records of the verified _absence_ of a species, the maps generated from SDMs should be interpreted as the relative environmental suitability for a species in a given location compared to other sites, rather than a probability of occurrence at a particular site --- though there are methods that allow for this translation @Smith2025LinkingRelative. 

In disease ecology, SDMs are often used to understand the current (and future) distribution of wildlife hosts and vectors (e.g. #cite(<Lippi2023TrendsMosquito>, form: "prose"), #cite(<Kopsco2022ScopingReview>, form: "prose"), #cite(<Bussieres-Fournel2025ClimateChange>, form: "prose")) and determine which areas have the potential to carry spillover risk. While SDMs can help in planning and management decisions (e.g., sampling prioritization), their interpretation can prove challenging for end users @Guisan2013PredictingSpecies. For example, predictions for each species (e.g. in a multi-host system) are frequently stacked (i.e. the suitability scores for each species are summed at each location) to understand their overlapping distribution. However, SDMs do not typically incorporate biotic interactions (for discussion of joint species distribution models, which incorporate co-occurrence data and are generally used to infer relationships rather than create spatially-explicit predictions, see #cite(<Pollock2014UnderstandingCooccurrence>, form: "prose") & #cite(<Wilkinson2021DefiningEvaluating>, form: "prose")), which means they cannot predict species abundance or the definitive interactions between species within a pixel where species are predicted to co-occur @Blanchet2020CooccurrenceNot. Without the influence of other species and other latent variables, SDMs likely overestimate the “true” suitable area (false positives) and overall number of co-occurring species. Additionally, SDMs typically use predictors that are treated as functionally static in time (e.g. WorldClim bioclimatic variables), which bakes in the assumption that species are in equilibrium with their environment @Milanesi2020IntegratingDynamic. In reality, species distribution patterns likely shift through time, varying with seasonal patterns, ecosystem disturbance, and resource availability @Milanesi2020IntegratingDynamic. Thus, when incorporating SDMs into a decision-making process, users must consider the plausibility of predictions and inherent uncertainty in the context of their particular expertise.

Further complicating the interpretation of SDMs are the many decisions made during data selection, pseudo-absence generation, model specification, and evaluation that proliferate through the outputs and generate errors @Merow2013PracticalGuide @Barry2006ErrorUncertainty. Even a poorly designed model can successfully be trained and generate predictions that seem plausible, but are ultimately as unreliable as the underlying model @Fourcade2018PaintingsPredict. Thus, there’s no replacement for collaboration between modelers with technical expertise and local decision-makers, with a clearly defined purpose and standardized reporting and documentation @Guisan2013PredictingSpecies @Araujo2019StandardsDistribution.

] <box1> 
