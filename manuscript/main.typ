#import "template.typ": *
#import "@preview/wordometer:0.1.4": word-count, total-words
#import "@preview/colorful-boxes:1.4.3": *
#import "authors.typ": *

#let title = "A protocol for biodiversity-informed wildlife disease surveillance"
#let abstract = "Land use and climate change threaten to increase the risk of spillover of zoonotic disease into human populations, creating a severe risk for public health. However, we typically lack actionable information about the prevalence of pathogens in wildlife populations for most of the globe. Even when this data exists, it is often sampled opportunistically and without guidance based on known presence of hosts of zoonotic pathogens. We present a protocol for integrating increasingly available data on host biodiversity into sampling priority for wildlife disease surveillance by enriching host species distribution models with pathogen prevalence data (if available). This protocol has the flexibility to adapt to different levels of data availability, but still makes sampling recommendations based on a principled understanding of host distribution and pathogen biology. We view this framework as a basis for integrating long-term biosurveillance and biodiversity monitoring programs."

#set par.line(numbering: "1")


#showtitle(title)
#showauthors(authors, institutions, corresponding_email)
#showabstract(abstract)
#show: word-count.with(exclude: (figure, bibliography))
#show figure.caption: set align(left)
#show figure: set block(breakable: true)

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

#show heading: it => {
  text(it.body, font: "Helvetica Neue", weight: 500, fill: rgb("#2e5385"))
  v(0em)
}
#show heading.where(level: 2): it =>{
    v(0.3em)
    text(it.body, size:1em, font: "Helvetica Neue", style: "italic", weight: 400, fill: rgb("#1d8265")) 
    v(0em)
}
#show heading.where(level: 3): it =>{
    v(0.3em)
    text(it.body, size:1em, font: "Helvetica Neue", style: "italic", fill: rgb("#444"), weight: 400) 
    v(0em)
}
#show cite: set text(fill: rgb("#528fd1"))


= Introduction 

Despite the well-established links between biodiversity and the potential for zoonotic pathogen spillover, wildlife disease surveillance programs rarely integrate biodiversity data to guide sampling efforts. The relationship between biodiversity and disease risk is complex, as both biodiversity itself and biodiversity loss can influence the presence of pathogens in context-dependent ways @Keesing2010ImpactsBiodiversity @Carlson2025a. Host and pathogen diversity are entwined through ecological and evolutionary processes @Plowright2017 @Carlson2025a, which together shape zoonotic spillover risk by affecting the prevalence of pathogens among wildlife that humans may encounter and from which infections may arise. Climate and anthropogenic land use change exacerbate this risk by modifying the ranges of hosts of zoonotic pathogens, leading to increased potential for cross-species transmission @Carlson2022, including spillover into human populations.

As a result, monitoring infectious diseases in wildlife is an important goal of several international agreements, which all fall under the One Health approach. The interlinkages between biodiversity and health are recognized in the Global Action Plan on Biodiversity and Health from the Convention on Biological Diversity (CBD) and the Pandemic Agreement of the World Health Organization (WHO). Specifically, the CBD encourages Parties to “[reinforce] planning and surveillance of biodiversity, including for wildlife habitats and zoonotic pathogen spillover risk, to better assess and address health and disease risks in order to manage wild species sustainably” (CBD/COP/DEC/16/19, Section IIIB). To meet this goal, countries must design monitoring systems that efficiently allocate resources to maximize their knowledge of zoonotic disease prevalence and spillover risk. Improving spillover risk estimation therefore requires a multifaceted approach that incorporates biodiversity surveillance with data on landscape change and host ecology --- including shifts in reservoir species’ habitats, community composition, and geographic distributions @Bell2025. 

The field of biodiversity monitoring has established many practices and tools that can directly inform biosurveillance efforts @Poisot2025b. One such core concept is the Biodiversity Observation Network (or BON): a set of monitoring locations designed to best capture the status and trends of biodiversity @Scholes2012. Ideally, BONs establish a set of spatial locations at which biodiversity data are collected in standardized formats, which can then be easily aggregated for the detection and attribution of biodiversity change @Gonzalez2023FrameworkDetection and be used to inform decision-making at local, regional, and ultimately global scales @Gonzalez2023GlobalBiodiversity. Here, we explore how the BON perspective would be useful to prioritize locations for which wildlife disease surveillance could be maximally informative, particularly when disease data are scarce. 

Given the cost of effective biosurveillance, strategic allocation of resources toward where, when, what, and how much to sample is imperative to ensure sampling effort yields new and useful information. Efforts have been made to think about biosurveillance sampling prioritization from both a statistical @Farver1985 @Nusser2008 and spatially explicit perspective @Andrade-Pacheco2020 @dumelle2022comparison. Here, we develop a context-dependent protocol for adaptive sampling to maximize the reduction of our uncertainty of the status of a given pathogen or host. The type of data (i.e. disease prevalence among a population of susceptible hosts, or occurrence records of host species) useful for a particular context depends on existing data availability and what drives our uncertainty of the hazard level of zoonotic disease prevalence. To optimally allocate sampling effort, both for choosing sampling locations and the type of data sampled, it is necessary to prioritize resource allocation towards both improving estimates of host presence and pathogen prevalence within host species. However, the locations for efficient sampling for these two goals are intrinsically different: to improve host range prediction we should go where we are most unsure about the host’s presence, whereas to improve prevalence we should  rather sample where we are already sure hosts are. 

Here, we establish a set of guidelines to use biodiversity data to prioritize sampling locations for wildlife disease surveillance. Notably, we select these sampling locations by accounting both for the distribution of biodiversity and current uncertainty about this distribution, as well as utilizing existing prevalence data if available. As applying this framework requires numerous decisions, we provide guidance about where it can be adapted to reflect local priorities, sampling constraints, and the availability of existing data on species presence and the prevalence of pathogens among host species. We present case studies illustrative of this framework in both data-poor (no prevalence data) and data-rich (spatially explicit estimates of prevalence) contexts on _Hanta_- and _Arenaviradae_ in rodents.

= The Role of Biodiversity Data in Biosurveillance 

Conceptualising zoonotic spillover requires distinguishing _hazard_ --- the active circulation of a pathogen in a host reservoir in a state conducive to transmission to humans --- from _risk_, which is the probability of that hazard being realised via spillover @Plowright2017. The intensity of this hazard fluctuates with shifts in host distribution, infection prevalence, and pathogen biology, all independent of human presence, and thus the proximity of reservoir hosts to human populations does not alone constitute a significant risk.

Spillover risk is actualized when a susceptible human is exposed to a zoonotic disease hazard, through processes governed by sociological, behavioural, and economic drivers at the human-wildlife interface @Friant2025. For instance, domestic or agricultural activities in environments contaminated by animal excreta elevate risk, whereas mitigations like improved sanitation or hygiene can reduce it towards zero. Monitoring the underlying hazard is critical because the conditions for exposure are not static. Future anthropogenic changes, such as urban encroachment or ecological shifts, like climate-driven species migration, can rapidly intensify human-wildlife contact, transforming a latent hazard into an immediate risk @Daszak2000 @Carlson2022.

It should be noted that zoonotic hazards cannot typically be directly estimated --- many unknown pathogens are expected to circulate in wild animal populations, which could be a source of novel human-infecting viruses. Still, the uncertainty created by this viral “dark matter” can be integrated into this protocol by adjusting the sampling priority for host species, using models to predict the most likely hosts of diseases that could jump into humans prior to an outbreak @Becker2022 @Cummings2025ViralEpidemic.

== When is biodiversity data relevant to guiding disease surveillance?

One of the leading ecological determinants of spillover risk is the density of infected hosts in a given location across all possible host species @Plowright2017. Unfortunately, this itself is often not actionable because it relies on spatially explicit data about abundance of --- and pathogen prevalence in --- all relevant hosts, which is rarely (if ever) available. Therefore, estimates of spillover risk are limited both by lack of knowledge on disease dynamics and the distribution and abundance of host species. 

To circumvent this common data limitation, we rely on a proxy that uses weighted estimates of local host species composition to prioritize locations for sampling and identify the type of sampling to conduct. While estimates based only on the predicted presence of host species do not fully capture the risk of disease transmission, they hint at the potential for pathogen presence, especially if several hosts can be weighted for e.g. pathogen competence, importance for transmission cycles, and the confidence a given species can effectively serve as a host. 

By contrast with a situation where all species are assumed to be equally important, prior knowledge about the relevance of different hosts for disease transmission can be incorporated in these models as a relative weight of host sampling priority. This data can come from phylogenetic similarity to known hosts, existing records of infection, or expert knowledge on the system in question @tseng2025viral. In the extreme case where there is no prior information to establish weights for host species, this approach lends itself to an adaptation of the well-established practice of estimating local host diversity by stacking predicted species distributions @Thuiller2015, which is particularly appropriate when using quantitative predictions of species presence @Grenie2020. Although stacking species distribution models (SDMs) is often biased towards higher estimates of local species richness @Calabrese2014, in our context this bias is unobjectionable as it leads towards more cautious recommendations that in the worst case overestimate the ranges of reservoirs, and the stacking of species-specific SDMs was shown to produce better estimates of richness compared to joint (multi-species) SDMs @Zurell2020. For a full set of guidelines for interpreting the outputs of SDMs, see #link(<box1>)[Box 1].

#pagebreak()
#colorbox(
  title: text(font: "Helvetica Neue", weight: 500)[Box 1: Guidelines for the interpretation of Species Distribution Models], 
  color: (
    fill: rgb("#f4fcfa"),
    stroke: rgb("#418b7f"),
    title: rgb("#101010")
  ),
  radius: 2pt,
)[  
  
Species distribution models (SDMs) are widely used tools in ecology that relate species occurrences with environmental variables to identify how species are geographically constrained by different ecological factors @Elith2011. SDMs can then be used to estimate the realized ecological niche of a species --- the environmental space (rather than geographic space) that is suitable for the species to occupy. Because occurrence data is often presence-only, meaning there are few/no records of the verified _absence_ of a species, the maps generated from SDMs should be interpreted as the relative environmental suitability for a species in a given location compared to other sites, rather than a probability of occurrence at a particular site --- though there are methods that allow for this translation @Smith2025. 

In disease ecology, SDMs are often used to understand the current (and future) distribution of wildlife hosts and vectors (e.g. #cite(<Lippi2023>, form: "prose"), #cite(<kopsco2022scoping>, form: "prose"), #cite(<bussieres-fournel2025climate>, form: "prose")) and determine which areas have the potential to carry spillover risk, although SDMs are not suitable to be applied directly to pathogens (#cite(<Carlson2020SpeciesDistribution>, form: "prose"), #cite(<Chipperfield2020>, form: "prose")). While SDMs can help in planning and management decisions (e.g., sampling prioritization), their interpretation can prove challenging for end users @Guisan2013. For example, predictions for each species (e.g. in a multi-host system) are frequently stacked (i.e. the suitability scores for each species are summed at each location) to understand their overlapping distribution. However, SDMs do not typically incorporate biotic interactions (for discussion of joint species distribution models, which incorporate co-occurrence data, see #cite(<Pollock2014UnderstandingCooccurrence>, form: "prose") & #cite(<Wilkinson2021DefiningEvaluating>, form: "prose")), which means they cannot predict species abundance or the definitive interactions between species within a pixel where species are predicted to co-occur @Blanchet2020. Without the influence of other species and other latent variables, SDMs likely overestimate the “true” suitable area (false positives) and overall number of co-occurring species. Additionally, SDMs typically utilize relatively coarse-grained temporal data that are functionally static in time, which bakes in the assumption that species are in equilibrium with their environment @Milanesi2020. In reality, species distribution patterns likely shift through time, varying with seasonal patterns, ecosystem disturbance, and resource availability @Milanesi2020. Thus, when incorporating SDMs into a decision-making process, users must consider the plausibility of predictions and inherent uncertainty in the context of their particular expertise.

Further complicating the interpretation of SDMs are the many decisions made during data selection, pseudo-absence generation, model specification, and evaluation that proliferate through the outputs and generate errors @Merow2013 @Barry2006. Even a poorly designed model can successfully be trained and generate predictions that seem plausible, but are ultimately as unreliable as the underlying model @Fourcade2018. Thus, there’s no replacement for collaboration between modelers with technical expertise and local decision-makers, with a clearly defined purpose and standardized reporting and documentation @Guisan2013 @Araujo2019StandardsDistribution.

] <box1> 


== What are host distribution models telling us about disease?

Directly applying SDMs to pathogens (such as viruses or obligate intracellular bacterial infections) can lead to results that don’t reflect reality @carlson2022global @Chipperfield2020. The circulation of these pathogens through an animal population can vary significantly over time and is influenced by host species behavior, immunology, and environmental factors, among other drivers. Over long time periods, pathogens are thought to circulate stably within their primary host populations throughout a landscape according to the contacts of individual members of a population. These pathogens can circulate through other susceptible animal populations if the two competent animal populations overlap in the landscape and the pathogen is shared. 

In places where no susceptible or potentially susceptible hosts are present, it is not possible to find infected individuals. Thus, the extent of a pathogen can be thought of as the aggregation of ranges for each host species, and a rough estimate of pathogen suitability can be obtained based on weighting host susceptibility to infection. The ensemble distribution then represents the landscape that the pathogen faces, which is not directly related to the physical environment, but is instead determined by the distribution of hosts in a landscape. 

= A Protocol for Biodiversity-informed Biosurveillance

Here we present a protocol for adaptive sampling of wildlife disease that adjusts sampling priority, both across space but also in terms of the type of data to collect (host occurrence vs disease prevalence), based on the availability of existing data. This provides a standardized framework for practitioners to decide where to sample, reevaluate priorities after data is collected, and design further samples to most efficiently inform the predicted status of zoonotic disease in wildlife.  

#figure(
  image("../plots/concept.svg"),
  caption: [*(a)*: The conceptual framework for optimizing sampling priority using host species distribution models, and pathogen prevalence data if available (left). *(b)*: The type of sampling that is most informative (pathogen prevalence vs. host occurrence), depending on a location's level of biodiversity dose (weighted host habitat suitability) and uncertainty about that dose. 
]
) <concept>
\
\
== Biodiversity dose, host weighting, and uncertainty

Throughout this manuscript, we refer to the estimate of local host composition that reflects both host richness and _a priori_ estimates of host importance as the “biodiversity dose” --- the ecological capacity for a location to maintain the pathogen (not a measure of the viral load in the environment). This captures both the likely presence of host reservoirs and their relative importance in disease transmission. Biodiversity dose is defined for a group of hosts, ideally encompassing the full set of hosts for a given pathogen, but also accounting for the subset of hosts of interest for a given monitoring program. 

In its simplest expression, the biodiversity dose is calculated as a _weighted_ average of SDM habitat suitability scores, but this still leaves the task of assigning species weights. If no information is available about species importance, a naive assumption of equal weights can be applied. If the capacity of a species to serve as a reservoir for a given pathogen is predicted with a model @Becker2022, the relative confidence in the prediction can be used as a weight. Hosts for a particular pathogen identified through methods with higher false-negative rates (e.g. serology vs. PCR) can be given lower weights. When there is uncertainty about the capacity of a host to be infected by a pathogen, methods that assign weights based on phylogenetic proximity to known hosts @Elmasri2020 can also provide usable information. Ideally, all aspects of the competence and relevance of hosts for disease transmission can be aggregated by experts to establish species weighting — we discuss the use of sensitivity analyses to better understand the consequences of these choices on the recommended configuration of the monitoring network in the final subsection of the protocol. 

Another crucial element for guiding sampling is the uncertainty about the distribution of the biodiversity dose, i.e. the weighted average of the uncertainty associated with the SDM for each host. There are several different ways to derive uncertainty from an SDM, and they tend to reflect different aspects of model fit --- see Box 3 for more information on the various forms of SDM uncertainty that can be used. High uncertainty points to an increased priority for sampling host occurrence to better refine the prediction of host ranges and overall biodiversity dose.

The information about biodiversity dose, dose uncertainty, and prevalence data (if available) can be used to make a map of sampling priority to decide where and what to sample --- see @concept\(b). It is critical to emphasize that while this diagram in Figure 1 is useful for planning the type of data to sample when prioritizing data collection, it cannot be used to infer risk or guide outbreak response. Instead, this protocol helps identify what type of sampling is most appropriate at a given stage of surveillance, based on host biodiversity distributions and assumptions about host importance. When the sampling objective is population surveillance or estimating pathogen prevalence within hosts, areas of high biodiversity dose (@concept\(b), Quadrants A & B) should be prioritized. However, when the goal is to improve host occurrence data, sampling should focus on areas of high dose uncertainty (@concept\(b), Quadrants B & D). 

Based on the goal of sampling, we can combine dose and uncertainty predictions into an overall priority score $P_x$ at each location $x$ as
 
$
  P_x = sum_i w_i (w_d s_(i x) + w_u u_(i x))
$

where $s_(i x)$ is the SDM score of species i at location $x$, $u_(i x)$ is the uncertainty of this score, $w_d$ is the weight for the biodiversity dose, $w_u$ is the weight for the uncertainty such that $w_u + w_d = 1$, and $w_i$ is the weight for species i, such that  $sum_i w_i = 1$.
\
\
== Integration of prevalence data into the monitoring network design

When available, data on pathogen prevalence can refine sampling priorities, but the best way to use this information for sampling is contingent on the spatial and temporal scale of the prevalence data and the goals of sampling. For example, continued surveillance of populations known to have high prevalence may be useful if the goal is to determine if this population serves as a persistent reservoir for a given pathogen, or if it was undergoing a transient outbreak when the existing data was collected. 

Alternatively, if existing prevalence data has high enough spatial coverage that it enables spatially explicit estimates of prevalence, sampling can most effectively improve spatial predictions of prevalence by targeting regions where the degree of prevalence is uncertain, but for which we are confident there is a high biodiversity dose (as we do in the second case study). This is useful for both precisely locating existing reservoirs across space to better map zoonotic hazard and spillover risk. For a full overview of potential methods to use depending on the scale of existing prevalence data see Supplemental Table 1, and for guidelines on how to target different forms of sampling regions based on their predicted dose, dose uncertainty, and predicted prevalence in Supplemental Table 2.

Another important consideration is the variability in the quality of prevalence data depending on the method used to detect infection. For example, serological diagnostic tests are less sensitive and specific than PCR, and in some cases are only available at a coarser taxonomic resolution. Further, fully representative sampling of wildlife populations is rarely feasible, so estimates of prevalence are obtained from a small subsample of the entire population of unknown size. Both of these factors induce error in the resulting prevalence estimates; however, there are numerous methods from population ecology for estimating abundance with imperfect detection (e.g. with N-Mixture Models #cite(<royle2004nmixture>, form: "prose")) and to infer population trends over time @kery2009trend. These can naturally be extended to the context of imperfect detection of both host abundance and pathogens within hosts in infectious disease modelling @direnzo2019diseasestructured.

== Sampling Point Selection

Once we have a sampling priority map, we want to use it to guide sampling site selection. There are many options for point selection algorithms rooted in sampling theory, which itself has a well-developed theory of spatial sampling. A full review of spatial sampling theory is outside of the scope of this paper, but see #cite(<Wang2012ReviewSpatial>, form: "prose") and #cite(<dumelle2022comparison>, form: "prose"), and for a comprehensive review of the general theory of sampling, see #cite(<Thompson2012Sampling>, form:"prose").

A crucial point here is that because we are interested in targeting areas of high priority, we are constrained to a set of point-selection methods for unequal probability sampling, meaning each location in space can have a distinct probability of being included in the sample. A naive approach would be to draw samples with probability directly proportional to priority value, but due to the autocorrelation associated with biodiversity dose and uncertainty, this would result in many sampling points clumped close together, which would likely provide redundant information, and be an inefficient use of sampling effort. This is a well-established issue in spatial sampling, and as a result there are many algorithms for selecting spatially balanced samples --- for example: Generalized Random Tessellation Sampling (GRTS; #cite(<stevens2004spatially>, form: "prose")), the Pivotal Method @grafstrom2012spatially, and Balanced-Acceptance Sampling (BAS; #cite(<robertson2013bas>, form: "prose")). These methods are all included in the BiodiversityObservationNetworks.jl package in Julia (which is maintained by the first author), which we utilize in the case studies.  

Recent results suggest that most site selection algorithms achieve equivalent performance for the same problem @Norman2025a, which allows users to pick an algorithm based on specific features, such as the support for auxiliary environmental data. For our purposes, we find BAS the most flexible method because it is both the most computationally efficient and achieves better spatial balance than the alternatives @robertson2013bas.

In practice, algorithms for generating spatially balanced samples targeted toward regions of high priority may not be skewed “enough” toward regions of high priority for a given use case, so we also tend to tilt the distribution of priority scores to make the values of high priority even higher. This is done using a function $T(P_x , alpha )$, defined as 

$
  T(P_x, alpha) = exp(alpha P_x)/ (1 + exp (P_x) )
$

where $P_x$ is the original priority at a given location $x$, and the parameter $alpha in (0, infinity)$ controls how much to tilt the adjusted map toward areas of high priority.

Another common method for designing samples is stratification. In the context of spatial sampling, stratification consists of splitting the spatial domain into discrete regions (corresponding with specific strata), which each contains a user-determined number of samples. In our case, stratification can be used to divide space into the different data collection regimes (sampling for host presences, prevalence or both; @concept\(b)) and distribute sampling effort across these strata based on the goals of the monitoring program, with spatially balanced point selection methods being used within each stratum. For further guidelines for interpretation and practical considerations when using selected sampling points, see *Box 2*.

#pagebreak()
#colorbox(
  title: text(font: "Helvetica Neue", weight: 500)[Box 2: Guidelines for the interpretation of recommended sampling locations], 
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
  @Altizer2004 @Cosgrove2008. Factors like migration and birth rates can affect these within-year cycles @vanDijk2014. Between-year fluctuations, particularly in rodents, are both dramatic and can also counter-intuitively inverse effect on populations’ prevalence @Davis2005. Here, we do not consider temporal variation in host abundances in our sampling recommendation protocol, since to do so would require orders of magnitude more data than we consider available, and would likely require the use of a different class of models/approaches. 

  #figure(
  image("../plots/annotated_sites.png"),
  caption: [Three selected sites for the South Korea case study (colored markers, left) overlayed on a map of South Korea with the 10 most populous cities shown in dark grey dots, and the major highway network shown in light grey (to get a sense of the accessibility of these sites). The images on the right correspond to the satellite imagery and nearby roads associated with each location, derived from OpenStreetMap, reflecting the heterogeneity within these pixels.]
) <box-figure>

] <box2>





== Sensitivity/Robustness Analysis

This protocol relies on several choices of weights based on the goals of sampling and _a priori_ estimates of host relevance for the pathogen(s) of interest. Because these values have a strong impact on sampling priority across space, assessing how sensitive the overall priority map is important to ensure the resulting priority is robust to small tweaks to the weights that do not reflect meaningful differences in host relevance. 

For a given set of weights, $bold(w)$, a natural way to assess the sensitivity of the overall priority is to examine the mean absolute change in the priority map, $P$, when a small random perturbation is applied to the weights. This gives us an interpretable metric of how sensitive the resulting priority map is at a given chosen weight value. We create a “nudged” version of the weights, $bold(w)’ = bold(w) + epsilon$, by adding noise $epsilon$, which is a vector of i.i.d. $"Normal"(0, sigma_("noise"))$ of the same length as $bold(w)$, and then renormalized (so $bold(w)'$ also sums to 1). We then construct the priority map, $P$, for the original weights, and the priority map $P'$ from the “nudged” weights. The overall sensitivity, $S$, of the priority map at a given weight value $bold(w)$ can be computed as the sum of the absolute difference across all locations x:  

$
  S(bold(w)) = sum_x |P_x - P'_x|
$

As an example, in Figure 2 we apply random noise to each possible value of the weights for the three host species we assess in the second case study. The ternary plot in the left of Figure 3 represents each possible combination of weights for three species. The color of each point is the sensitivity $S(bold(w))$ for noise drawn from standard Normal distribution with $sigma_("noise")=0.01$ across 100 samples. 


If the selected weights for a given use case are near the "warm" values on the heatmap, end users should consider the additional caveat that the selected weights are very sensitive to small changes, and ensure the selected weight values are   reflective of meaningful biological differences in host relevance.

#figure(
  image("../plots/sensitivity.png"),
  caption: [The sensitivity of the overall priority map to adjusted weights. The color of each point on the left ternary diagram is proportional to the sensitivity of each point to nudging weights (see main text). Each panel on the right corresponds to the overall priority map at weight values in the corresponding color on the left.]
) <sensitivity>



#pagebreak()
#colorbox(
  title: text(font: "Helvetica Neue", weight: 500)[Box 3: Guidelines for the refinement of various methodological steps], 
  color: (
    fill: rgb("#fcf7ff9d"),
    stroke: rgb("#866999"),
    title: rgb("#101010")
  ),
  radius: 2pt,
)[ 

The workflow described in this manuscript involves many steps, and most of them call for decisions made by end users. In this box, we go through the different steps, and highlight key methodological considerations.

*List of hosts*: the list of hosts can be assembled from biodiversity data (IUCN range maps, in-country checklist, GBIF data), or from past wildlife disease data (list of competent hosts from e.g. VIRION @carlson2022global, literature surveys, or previous surveillance programs). The list of hosts may be arbitrarily filtered or expanded to reflect local priorities.

*Hosts weights*: weights of hosts used in the biodiversity dose calculation should essentially serve as a quantification of their expected contribution to disease transmission, and will almost always be a multi-factorial compromise between e.g. strength of evidence (serology is less credible than PCR), and competence or transmission risk (potentially estimated using phylogenetic similarity to well-known hosts). Weights can also ultimately attempt to capture the risk of transmission to human populations (e.g. by weighting synanthropic species higher).

*Host SDMs*: the usual recommendations about the training and validation for SDMs apply throughout this pipeline (see #cite(<zurell2020standard>, form: "prose")), and in particular the choice of predictor data, spatial extent, the quality control of occurrences, and the selection of predictive variables are important.

*Uncertainty*: There are many ways to quantify SDM uncertainty. Some models (like Generalized Linear Models and the form of a Boosted Regression Tree we use in the case studies) have estimates of variance built into the model structure. Absent this, a common alternative is measuring the variance of the predicted suitability score at each site across many cross-validation folds, or methods using conformal prediction @Poisot2024ConformalPrediction. An alternative is using Bayesian Additive Regression Trees (BART; #cite(<Carlson2020EmbarcaderoSpecies>, form:"prose")), which directly obtains samples from parameters, and which can thereby be aggregated into uncertainty in predicted suitability. 

*Prevalence*: Relevant factors for choosing what prevalence data to incorporate are the minimum number of individuals sampled, the time since data was collected, and the form of test used. These can each impact the reliability of the prevalence estimate, and therefore utility of prevalence data. 

*BON design*: Many algorithms exist for spatially balanced sampling with unequal inclusion probabilities, and the ability of these algorithms to handle various forms of auxiliary data should be considered (see #cite(<Norman2025a>, form:"prose")). It is necessary to make the typical considerations when planning sampling: the accessibility of a given site, the cost-effectiveness of sampling given how much time is required to reach a location, the ability to access private land, the sovereignty of indigenous land --- these can all be used to further adjust the inclusion probability. 

*Sensitivity analysis*: Assessing the sensitivity of the overall priority map to both selected species weights (as discussed in the section on sensitivity), as well as the weights toward uncertainty and dose are all key factors to ensure the priority map is robust to small changes in weighting. Further checks on the response of the overall priority map to the inclusion/exclusion of host species (particularly those for which there is no direct evidence they can host a particular pathogen) should also be considered.

] <box3>
#pagebreak()

= Case Studies 

To illustrate this approach, we consider two case studies on _Arena_- and _Hantaviridae_ and their rodent hosts, using the global database compiled by @simons2025protocol. The first case study is on rodent hosts of _Arenaviridae_ in India, and reflects the scenario where there is no prevalence data available. The second is for South Korea, where we have sufficient prevalence data on _Hantaviridae_ to make spatially explicit predictions of prevalence and use this to guide sampling. 

== Creating Dose and Uncertainty Maps

For both case studies, we use the same SDM methodology outlined here. We emphasize that the specific approach to building SDMs is not the focus --- sampling prioritization is only going to be as good as the SDM, and we present guidelines for the interpretation of SDM outputs in #link(<box1>, "Box 1"). 

We obtain occurrence records for each rodent host species from the Global Biodiversity Information Facility (GBIF). As predictor variables, we use the 19 bioclimatic variables from CHELSA @karger2017climatologies at 1km x 1km resolution for South Korea, and from WorldClim @fick2017worldclim at 2.5 arcminute resolution (approximately 5km x 5km at the equator) for India (to keep the raster size reasonable). We generate pseudoabsences @Barbet-Massin2012 using background thickening @Vollering2019, with a buffer around each occurrence record where no pseudoabsences can go (25km for India, 8km for South Korea). 

Models were trained using SpeciesDistributionToolkit.jl @Poisot2025e and EvoTrees.jl @Desgagne-Bouchard2025EvovestEvotreesjl in Julia. Specifically, we use a Boosted Regression Tree (BRT; #cite(<Elith2008WorkingGuide>, form: "prose")) with a Gaussian loss metric, meaning the value of each node in the tree is fit to a Gaussian using maximum-likelihood estimation. Therefore, the BRT provides both a suitability prediction score, and an uncertainty value associated with it, for each pixel.

In order to standardize the dose values for each species, we compute the empirical cumulative distribution function (ECDF) on raw output predictions, ensuring all prediction layers are on the same scale. 

== Case Study 1: Designing a network without prevalence data

Our first case study considers the situation where we are primarily limited by a lack of prevalence data. This case study was inspired by #cite(<Rodrigues1978SerologicalSurvey>, form: "prose"), published 9 years after the Lassa virus (_Mammarenavirus lassaense_) was first identified in Nigeria. In the immediate aftermath of the discovery of a virus of severe public health importance that circulates in hosts with widespread distributions, #cite(<Rodrigues1978SerologicalSurvey>, form: "prose")  were the (self-identified) first to sample for Lassa virus outside of Africa. At the time, Lassa itself was difficult to distinguish in clinical settings from other infectious diseases present in India (e.g. typhoid fever, malaria), and therefore it was not unreasonable to suspect it may be present in rodent populations in India.    

Today we know that Lassa is confined to Western Africa, but this case study serves as an example for how we use biodiversity data to prioritize sampling in the context of a recently identified virus without existing reliable information about its geographic extent, prevalence, or known reservoir species. We construct a priority map for sampling using the host rodents tested by #cite(<Rodrigues1978SerologicalSurvey>, form: "prose") : _Funambulus tristriatus_, _Hystrix indica_, _Rattus rattus_, and _Suncus murinus_. To generate a priority map for sampling, we follow the SDM procedure outlined above, and for the sake of example assign species weights at random. 

@india\(a) shows a bivariate map of both the dose and uncertainty aggregated across host species. This highlights three large-scale regions of interest: the south, where the model is confident that many host species are present, the northwest, where the model is uncertain about many or all species, and the northeast, where the model is also relatively confident that host species are present. We say a region in the top 50% of uncertainty for a given species is for “discovery sampling”, and a region in the top 50% of predicted suitability is for “prevalence sampling”. In @india\(b) we see the total number of species that fall into each of these categories across space. @india\(c) shows our priority map, with both $w_u$ and $w_d$ equal to 0.5 (reflecting equal priority toward discovery and surveillance). The white points reflect sampling locations generated using Balanced Acceptance Sampling applied to the logistically tilted priority map with $alpha = 5$ (see Sampling Point Selection subsection above). The teal markers reflect the four historical sampling locations from #cite(<Rodrigues1978SerologicalSurvey>, form: "prose"). Finally, in @india\(d), we see the host species that contributed the most to the priority score across space (computed as the maximum value of $w_i (w_u U_x + w_d  P_x)$ at location $x$ across species $i$). 


#figure(
  image("../plots/india.png"),
  caption: [The case study for India. (a) Biodiversity dose + uncertainty bivariate plot. (b) Number of hosts in prevalence-regime (top 50% of habitat suitability for that species) vs. discovery (top 50% of uncertainty for that species), (c) Sampling priority + BON. X’s are historical sampling locations from #cite(<Rodrigues1978SerologicalSurvey>, form: "prose") (d) The largest host contributor to priority]
) <india>

== Case Study 2: Accounting for prevalence data in network design

In contrast to the prior example where very little is known about the prevalence of a given pathogen, we now demonstrate the protocol for the surveillance of two hantaviruses in South Korea (_Orthohantavirus hantanense_ & _Orthohantavirus seoulense_), which are both of RNA-virus order Elliovirales @Vaheri2013, and are two of the main causative agents of hemorrhagic fever with renal syndrome (HFRS; #cite(<Park2025>, form: "prose")). Typically, the hosts of these viruses are rodents and shrews, and in our study we focus on three small-mammal hosts: _Rattus norvegicus_, _Apodemus agrarius_, and _Crocidura lasiura_. These host-virus pairs were selected as the database contained prevalence records across the country at a broad enough spatial scale to enable spatially explicit prediction of prevalence. The data comes from the following studies: #cite(<gu2011genetic>, form: "prose"), #cite(<kim2017hantavirus>, form: "prose"), #cite(<kim2007seroepidemiological>, form: "prose"), #cite(<kim2020active>, form: "prose"), #cite(<klein2015hantaan>, form: "prose"), #cite(<lee1978isolation>, form: "prose"), #cite(<lee2017dynamic>, form: "prose"), #cite(<no2016detection>, form: "prose"), #cite(<park2022portable>, form: "prose"), #cite(<ryou2011prevalence>, form: "prose"), #cite(<sames2009ecology>, form: "prose"), #cite(<seo2022emerging>, form: "prose"), and #cite(<song2009characterization>, form: "prose").

Prevalence across the extent is estimated using Gaussian Process regression (also known as kriging) for each host-virus pair with at least 5 populations sampled across space. This is done using a Squared-Exponential covariance kernel with a fixed length scale parameter using GaussianProcesses.jl @fairbrother2022gaussianprocessesjl, and the hyperparameters for the model are optimized using the L-BFGS algorithm in Optim.jl @mogensen2018optim. This provides both spatial estimates of prevalence, but crucially, spatially explicit uncertainty associated with this prevalence estimate.


@korea shows both the bivariate plot of biodiversity dose and dose uncertainty (@korea\A) and dose vs. prevalence uncertainty (@korea\B). To isolate regions where we are confident hosts are present for prevalence sampling, and regions of high host presence uncertainty for discovery sampling, we stratify the spatial extent into discrete regions that prioritize prevalence sampling and occurrence sampling as follows: the region that is both in the top 50% of biodiversity dose and top 50% of prevalence uncertainty is the stratum for prevalence sampling, and the region in the top 50% of dose uncertainty is the stratum for discovery sampling (@korea\C). Spatial regions that fall into both categories are then targets for both forms of sampling. The choice of 50% quantiles is arbitrary for this example. The sampling priority within each of these strata can be seen in @korea\D. The overall priority map can be seen in @korea\E,  along with the points sampled using Balanced Acceptance Sampling with a tilting value of $alpha = 2.0$, and the type of sampling that should take place at each location. 

This case study demonstrates how we can integrate host SDMs with existing prevalence information to design sampling strategies and split geographic space into discrete regions where each type of sampling (occurrence and prevalence) should take place.



#figure(
  image("../plots/korea.png"),
  caption: [The case study for South Korea. (*A*): Bivariate plot of dose (increasing green) and dose uncertainty (increasing blue). (*B*): Bivariate plot of dose (increasing green) and prevalence uncertainty (increasing orange; measured as the sum of estimated variance from kriging across each host-virus pair). (*C*): The spatial strata for which occurrence, prevalence, and both forms of sampling should occur (see main text for how these are computed). (*D*) The priority scores for regions where discovery sampling is a priority (increasingly blue) and prevalence sampling is the priority (increasingly purple) in their respective strata. (*E*) The total priority map with sites selected using Balanced Acceptance Sampling. The marker shape of each site corresponds to the type of sampling that should occur there, based on the stratum in (C) that each point falls within.]
) <korea>



= Discussion 

Here, we have presented a protocol for adaptive sampling of wildlife disease based on weighted predictions of host richness and uncertainty, and guidelines for integrating prevalence data into sampling prioritization if available. The strength of this approach is prioritizing sampling locations for wildlife disease surveillance by accounting for varying availability of prevalence and host competence data, allowing it to be tailored to any particular stage of developing biosurveillance programs. When enough prevalence data are known to allow the spatially explicit modelling of prevalence @Gras2018, relying on biodiversity data becomes a lower priority. Still, because spatially explicit modelling of prevalence is mostly feasible at small spatial scales, and because it is important to keep track of larger-scale changes in species distribution @Carlson2022, we anticipate that the design of a biosurveillance network that is aligned with biodiversity data collection should become a standard practice. We discuss further methodological refinements that can spur future directions in #link(<box3>, "Box 3").

While this protocol relies on the spatial distribution of suitable habitat for pathogen reservoirs, our approach is not designed to capture temporal variance in virus prevalence or transmission dynamics, nor is it particularly suited for guiding time-series based optimization of surveillance efforts. Instead, this framework identifies areas with suitable environmental conditions for pathogen hosts, but does not predict when virus prevalence, or spillover risk, is increased within those areas. To fully maximize the effectiveness of sampling, users must integrate domain expertise regarding viral ecology, seasonal transmission patterns, and host-pathogen dynamics to interpret model outputs appropriately. Additionally, the static nature of SDMs may not adequately reflect rapidly changing environmental conditions, anthropogenic landscape modifications, or evolving host community structure that could influence reservoir distributions or pathogen prevalence in time and space. 

There are numerous future directions for incorporating more sophisticated forms of SDMs to help alleviate these challenges: temporally-explicit SDMs are possible but require much more data at a high temporal resolution, which are typically not available for most species @Anselmetto2025. Similarly, mechanistic SDMs, which directly integrate biological processes into habitat suitability predictions, require far more data than is typically available. Another limitation of conventional correlative SDMs is that habitat suitability does not directly translate to species occupancy. Species may not be present in suitable habitat due to competition or dispersal limitation, and this can impact viral transmission if a host species is replaced by a less competent host (e.g. Lassa fever is reported less in cities where _Mus musculus_ has displaced primary reservoir species). Joint SDMs (JSDMs), which model species composition across an entire species pool, thereby accounting for species interactions, are another possible frontier, but JSDMs require explicit co-occurrence data, typically not attainable from platforms like GBIF. Generally, these more informative SDMs are not possible without extensive, long-term monitoring, which is why the goal of this protocol is to get these long-term monitoring programs off the ground with efficient sampling in the low data context.

However, an increasing amount of pathogen prevalence data is only actionable when individual test outcomes, both positive and negative, are shared in a way that is both standardized and interoperable with the commonly used representations of biodiversity data @Schwantes2025. Further integration of open biosurveillance and biodiversity monitoring systems is essential to meet our common goals @Poisot2025b. Better biosurveillance itself is not spillover prevention, but it’s a necessary component of it to make prevention possible --- it’s the best we can do in the context of limited information on prevalence --- but a full approach to prevention of zoonotic spillover is broader than this alone @Ramakrishnan2023. This protocol represents a necessary first step toward widespread biosurveillance integrated with global biodiversity monitoring systems @Gonzalez2023GlobalBiodiversity to ensure a just planetary future. 

*Software and Data Availability*: The software and data used for this analysis is available in #link("https://github.com/gottacatchenall/OptimalSamplingBiosurveillance")[#text(fill:blue)[this]] GitHub repository.

*Acknowledgements*: MDC was funded through an IVADO postdoctoral fellowship. TP was funded by the Wellcome Trust (223764/Z/21/Z); TP was supported by the US National Science Foundation (grant no. DBI 2213854). AB was funded by the CANOPY CIHR training program. CBB was supported by an NSERC PGSD (no. 56749). HR was funded by an NSF GRFP (grant no. DGE-2139841). BK was funded by CIHR DRA (BCB 201224). DS was funded by a joint NSF‐NIH‐NIFA EEID Award \#2208034


#showbibliography("refs.bib") 

#pagebreak()
= Supplemental Material

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
    [In the ideal case, where there is both high spatial and temporal coverage, hierarchical Bayesian methods from the spatial context can be extended to include time-series components --- see @Zhou2008EwmaSmoothing @boulieri2020bayesian @lee2014cluster. Similarly, ecological models that incorporate imperfect detection (e.g. N-mixture models) can be extended to spatiotemporal context @zhao2017spatially, which can then be extended to include the hierarchical nature of pathogen detection @direnzo2019diseasestructured.], 


)] <prev-data-types>

#include "big_chungus_table.typ"

