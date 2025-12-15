#import "@preview/colorful-boxes:1.4.3": *

#set text(lang: "en", size:10pt)


#colorbox(
  title: text(font: "Helvetica Neue", weight: 500)[Box 2: Guidelines for the refinement of various methodological steps], 
  color: (
    fill: rgb("#fcf7ff9d"),
    stroke: rgb("#866999"),
    title: rgb("#101010")
  ),
  radius: 2pt,
)[ 

The workflow described in this manuscript involves many steps, and most of them call for decisions made by end users. In this box, we go through the different steps, and highlight key methodological considerations.

*List of hosts*: the list of potentially competent hosts can be assembled from biodiversity data (IUCN range maps, in-country checklist, GBIF data), or from past wildlife disease data (e.g. from VIRION @Carlson2022GlobalVirome, literature surveys, or previous surveillance programs). The list of hosts may be arbitrarily filtered or expanded to reflect local priorities.

*Host weights*: weights of hosts used in the biodiversity dose calculation should essentially serve as a quantification of their expected contribution to disease transmission, and will almost always be a compromise between strength of evidence (e.g. serology, PCR and pathogen isolation provide different levels of confidence in host species competence), and transmission risk (potentially estimated using phylogenetic similarity to well-known hosts). Weights can also ultimately attempt to capture the risk of transmission to human populations (e.g. by weighting synanthropic species higher).

*Host SDMs*: the usual recommendations about the training and validation for SDMs apply throughout this pipeline (see #cite(<Zurell2020StandardProtocol>, form: "prose")), and in particular the choice of predictor data, spatial extent, the quality control of occurrences, and the selection of predictive variables are important.

*Uncertainty*: There are many ways to quantify SDM uncertainty. Some models (like Generalized Linear Models and the form of a Boosted Regression Tree we use in the case studies) have estimates of variance built into the model structure. Absent this, a common alternative is measuring the variance of the predicted suitability score at each site across many cross-validation folds, or methods using conformal prediction @Poisot2024ConformalPrediction. An alternative is using Bayesian Additive Regression Trees (BART; #cite(<Carlson2020EmbarcaderoSpecies>, form:"prose")), which directly obtains samples from parameters, and which can thereby be aggregated into uncertainty in predicted suitability. 

*Prevalence*: Relevant factors for choosing what prevalence data to incorporate are the minimum number of individuals sampled, the time since data was collected, and the form of test used. These can each impact the reliability of the prevalence estimate, and therefore utility of prevalence data. 

*BON design*: Many algorithms exist for spatially balanced sampling with unequal inclusion probabilities, and the ability of these algorithms to handle various forms of auxiliary data should be considered (see #cite(<Norman2025SiteSelection>, form:"prose")). It is necessary to make the typical considerations when planning sampling: the accessibility of a given site, the cost-effectiveness of sampling given how much time is required to reach a location, the ability to access private land, the sovereignty of indigenous land --- these can all be used to further adjust the inclusion probability. 

*Sensitivity analysis*: Assessing the sensitivity of the overall priority map to both selected species weights (as discussed in the section on sensitivity), as well as the weights toward uncertainty and dose are all key factors to ensure the priority map is robust to small changes in weighting. Further checks on the response of the overall priority map to the inclusion/exclusion of host species (particularly those for which there is no direct evidence they can host a particular pathogen) should also be considered.

] <box2>