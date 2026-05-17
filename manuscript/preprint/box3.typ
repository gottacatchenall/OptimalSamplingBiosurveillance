#import "@preview/colorful-boxes:1.4.3": *

#set text(lang: "en", size:10pt)


#colorbox(
  title: text(font: "Helvetica Neue", weight: 500)[Box 3: Guidelines for the interpretation of recommended sampling locations], 
  color: (
    fill: rgb("#f5f5fc65"),
    stroke: rgb("#566a8d"),
    title: rgb("#101010")
  ),
  radius: 2pt,
)[ 
  We want to emphasize that sampling points, though given as exact (longitude, latitude) coordinates, are not to be interpreted without critical thought. Since the grain at which the underlying SDMs are operating is often significantly larger than the scale at which one may have to decide where to place a sampling trap, selected points should be thought of as most useful as the centroids of "local" regions (where "local" corresponds to the pixel size of the SDMs used for prioritization) wherein sampling effect could best improve knowledge about the ecological settings of viral transmission. To make this point explicit, if the model suggests sampling for _Mus musculus_ in an urban area and the coordinate itself is in the middle of the center of a highway, the model’s output should be interpreted as indicating the general area where sampling should be done (see below figure). 

  For example, @box-figure shows three selected sampling sites for the second case study on South Korea, with the satellite imagery corresponding to the cell in the priority map that was selected. As the original climate layers have pixels that are approximately 1 $"km"^2$, there is considerable heterogeneity within these pixels. Experts on particular taxa should assess the scope of the local region near each sampling point, and use it to guide the practical logistics of where effort (e.g. traps) are most likely to result in useful information.

  Beyond just the spatial scale at which recommendations are made, temporal variation in both prevalence and host occurrence is a well-known driver of disease dynamics, often driven by environmental variation tied to seasonality 
  @Altizer2004SeasonalDynamics @Cosgrove2008SeasonalVariation. Factors like migration and birth rates can affect these within-year cycles @vanDijk2014JuvenilesMigrants. Between-year fluctuations, particularly in rodents, are both dramatic and can also counter-intuitively inverse effect on populations’ prevalence @Davis2005FluctuatingRodent. Here, we do not consider temporal variation in host abundances in our sampling recommendation protocol, since to do so would require orders of magnitude more data than we consider available, and would likely require the use of a different class of models/approaches. 

  #figure(
  image("../../plots/annotated_sites.png", width: 310pt),
  caption: [Three selected sites for the second case study (colored markers, left) overlayed on a map of South Korea with the 10 most populous cities shown in dark grey dots, and the major highway network shown in light grey (to get a sense of the accessibility of these sites). Right: Satellite imagery, derived from OpenStreetMap, for each site.]
) <box-figure>

] <box3>
