#import "template.typ": *
#import "@preview/colorful-boxes:1.4.3": *

#show table.cell.where(y: 0): set text(size: 0.7em)
#show table.cell.where(x: 0): set text(weight: "medium")
#show table.cell.where(y: 0): set text(weight: "bold")
#show table.cell: it => {
  if it.y == 0 or it.x == 0 {
    set text(white, font: "IBM Plex Sans")    
    set align(center)
    strong(it)
  } 
  else {
    it
  }
}
//#show table.cell: set par(justify: false)
#show table.cell: set text(size: 0.75em,)
#show table.cell.where(x: 0): set text(weight: "medium")
#show table.cell.where(y: 0): set text(weight: "bold", size: 1.5em)
#show table.cell: it => {
  if it.y == 0 {
    set text(white)    
    it
  } 
  else {
    it
  }
}

= A protocol for biodiversity-informed wildlife disease surveillance -- Supplemental Material

== Example Sensitivity Analysis 

As an example of weight sensitivity analysis, in @sensitivity we apply random noise to each possible value of the weights for the three host species we assess in the second case study. The ternary plot in the left of @sensitivity represents each possible combination of weights for three species. The color of each point is the sensitivity $S(bold(w))$ for noise drawn from standard Normal distribution with $sigma_("noise")=0.01$ across 100 samples. 


#figure(
  image("../../plots/sensitivity.png", width: 110%),
  caption: [The sensitivity of the overall priority map to adjusted weights. The color of each point on the left ternary diagram is proportional to the sensitivity of each point to nudging weights (see main text). Each panel on the right corresponds to the overall priority map at weight values in the corresponding color on the left.]
) <sensitivity>

== Methods for incorporating prevalence data with varying levels of spatial and temporal coverage

#figure(numbering: "1", caption: [Methods for incorporating prevalence data with varying levels of spatial and temporal coverage])[
#table(
  columns: (12%, auto, auto),
  inset: 10pt,
  stroke: none,
  fill: (x, y) => 
    if y > 0 and x > 0 {
      if calc.even(y) { 
          gray.lighten(90%)
        } 
      } 
    else { 
      rgb("#2c7a5eff") 
    },  
  align: (x, y) => {
    if x > 0 and y > 0 {
      left 
    }
    else {
      horizon 
    }
  },
  table.header(
    [], [Low Spatial Coverage], [High Spatial Coverage],
  ),

  [Low Temporal Coverage], 
    [With little, but some prevalence data, it is difficult to use it to inform sampling priority, but because there is _some_ existing data, effort could be shifted toward prevalence sampling to supplement existing data. Whether this means returning to the same sites (to generate better time-series) or diversifying spatial locations (to improve spatial prediction) is dependent on the overall goals of sampling.], 
    [When there are many spatial locations with prevalence data, but each has few/one time point, the prevalence component of spatial priority can be obtained by using kriging (as we do in the second case study),  also known as Gaussian process regression. There is a rich literature on adaptive sampling specifically for improving estimates for kriging @Fuhg2021StateoftheartComparative, which could be incorporated into the prevalence component of priority.

    Alternatively, Bayesian hierarchical methods can be used if the goal is to quantify the effect of the auxiliary variables on the prevalence of disease. This can be combined with kriging as Gaussian Processes naturally can be converted to the Bayesian context.], 

  [High Temporal Coverage], 

    [Time-series analysis, e.g. to determine if a given sampled location is either a persistent reservoir or underwent a transient outbreak. As there is not enough spatial prevalence data to make spatially explicit predictions, the proximity to high or low prevalence areas can be used to weigh priority (see the final subsection of the protocol).],
    [In the ideal case, where there is both high spatial and temporal coverage, hierarchical Bayesian methods from the spatial context can be extended to include time-series components --- see @Zhou2008EwmaSmoothing @Boulieri2020BayesianMixture @Lee2014ClusterDetection. Similarly, ecological models that incorporate imperfect detection (e.g. N-mixture models) can be extended to spatiotemporal context @Zhao2017SpatiallyExplicit, which can then be extended to include the hierarchical nature of pathogen detection @DiRenzo2019DiseasestructuredNmixture.], 


)] <prev-data-types>

#include "big_chungus_table.typ"

#showbibliography("refs.bib") 
