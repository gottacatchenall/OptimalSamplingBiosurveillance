#set text(font: "Libertinus Sans")
#show math.equation: set text(font: "Libertinus Math")

#let response(ed: false, body) = block(fill: color.hsl(195deg, 30%, 98%), inset: 1em, radius: 0.05em, width: 100%, stroke: (left: 0.2mm + color.hsv(195deg, 100%, 45%)))[
  #if ed == true {
    text(weight: "bold")[Comments for the editor]
    linebreak()
    linebreak()
  } else {
    text(weight: "bold")[Response]
    linebreak()
  }
  #text(font: "Libertinus Serif")[#body]
]


// Editing marks
#let add(body) = text(fill: rgb(0, 100, 0))[#body]
#let change(body) = text(fill: rgb(0, 100, 100))[#underline(body, stroke: rgb(0, 90, 90))]
#let cut(body) = text(fill: rgb(150, 150, 150))[#strike(body, stroke: rgb(100, 0, 0))]
#show "TK": text(weight: "bold", fill: rgb("#e08619"))[TK]

#response(ed: true)[
  The changes are also presented in a track changed manuscript uploaded as a document for reviewers, which indicates #add[additions], #change[changes], and #cut[deletions] to the text resulting from comments made by both reviewers.
]


= Reviewer 1

I have all of my comments on the attached PDF. But to summarize the strengths of the manuscript:

- I appreciate the idea of identifying sampling priorities while explicitly incorporating spatial uncertainty. Overall, I see this as a promising step forward in environmental risk mapping.
- While species weighting always involves assumptions and caveats, I applaud the authors for trying to bring more nuance to weighing decisions (eg serology vs PCR). In many modeling frameworks, weighting decisions are left to the modeler’s discretion and can feel somewhat ad hoc. The effort to provide a more structured approach is commendable.
- I found the distinction between sampling location and type of data sampled to be particularly insightful
- The conceptual distinction between hazard and risk is also useful and clearly framed

#response[We thank the reviewer for their kind words summarizing the manuscript's strengths, and their thorough comments and suggestions in response to the manuscript, which we believe have strengthened the manuscript considerably. ]

I was able to clone the repository and run part of the code. Without troubleshooting or
actively fixing errors, I was able to run the India case study but not the South Korea case
study. I initially suspected the issue might relate to OSM object, but it appears to be
unrelated. However, I am primarily an R user, so this could reflect my own skillset rather
than an issue with the repository itself. If may be helpful for the authors to ask another
primarily R user to test the workflow on a different machine.

#response[We believe that the reviewer's issue is likely related to an update that occurred from Julia 1.11 (which is the version of Julia in which the original code was developed) and Julia 1.12. A minor change to the LinearAlgebra package caused a function used in kriging to break.

We have both addressed this issue and pinned explicit package versions in the Julia environment file (Package.toml) to improve reproducability in the future, and added explicit instructions in the README.md about versioning requirements to avoid incidents like this in the future. 
]

== Major Comments (General Feedback)

Overall, I found this article to be well written, visually appealing, and clearly motivated by policy relevance in the Introduction. 

However, I repeatedly found myself getting stuck on the framing that BON – and this protocol in particularly – is the primary or only framework addressing these questions. There has been substantial effort & investment in standardized biodiversity sampling, and there have also been many attempts to link biodiversity data with disease ecology indices (e.g., pathogen prevalence, pathogen load, host community structure). 

Because the manuscript does not sufficiently situate itself within the broader literature or acknowledge earlier related work (outside of the co-author list), the framing occasionally felt somewhat narrow. This was frustrating at time because in most other respects the manuscript is quite strong. 

I recognize that several authors on this paper have made important contribution to the biodiversity-disease space (and I enjoy these papers!). However, this is a large and active field with many research groups working on related problems. I would encourage the authors to broaden the citation base so that the manuscript better reflects that broader community. A brief scan of the reference list suggests a substantial number of citations include members of the author team (even with the “et al” truncation). 

#response[

We have restructured part of the introduction to introduce both existing biosurveillance and long-term biodiversity monitoring programs, and highlighted that the novel perspective introduced here is the potential to integrate these monitoring systems. 

We have also added numerous additional citations throughout the manuscript, but particularly in the introduction.

We are open to any further references that the reviewer believes we have missed. 

]

Similarly, the manuscript highlights BON as the primary monitoring framework, despite the existence of many other long-term biodiversity monitoring networks. While BON is an important example, it would be more appropriate to frame it as one of several efforts working toward standardized biodiversity monitoring. This work builds on substantial existing foundation, and acknowledging the broader literature would strengthen the manuscript.

#response[
  We have added a paragraph in the introduction that covers various long-term biodiversity monitoring projects (BBS, LTER, NEON, TERN), the role of GBIF and data standards like Darwin Core, the Essential Biodiversity Variable (EBV) and Essential Ecosystem Service Varaible (EESV) frameworks, and the GBF, and how these can support biosurveillance.
]

== Specific Comments

=== _Introduction_

==== Comment 1:

More consideration could be given to citation choices in the opening paragraph of the Introduction. While the cited research groups have contributed significantly to this field, it would be beneficial to diversify the research groups cited in the opening paragraph, given the breadth of work in biodiversity-disease research.

#response[We appreciate the comment and have added numerous references to the introduction to widen the citation-base. We are happy to include any additional citations that the reviewer believes are missing. 
]

==== Comment 2: 

After reading lines 49-57, I wondered who would realistically fund the surveillance and standardization systems being proposed. Are these cost estimates or funding mechanisms that could support implementation of the monitoring protocol suggested here? While these logistical considerations may be too detailed for the Introduction, I hope the Discussion includes a more substantial consideration of the trade- offs between ideal surveillance systems and what is practical to implement. Treaties and global initiatives may call for these systems, but the protocol proposed here implies substantial additional effort (at the individual and institutional levels), and it would strengthen the paper to acknowledge the logistical challenges users may face.

#response[We have added an additional paragraph in the discussion on funding and integration with a Global Biodiversity Observing System @Gonzalez2023GlobalBiodiversity and its funding model.]

==== Comment 3: 

The manuscript focuses heavily on BON, but there are many long-term
monitoring networks conducting similar work. It would be useful to frame this approach as
building on existing monitoring infrastructure rather than positioning BON as the primary
monitoring network. For example, within North America alone there are several well-
established networks including US LTER, US NEON, North American Breeding Bird Survey.
Highlighting similar networks would help acknowledge the substantial investments that
have already been made in biodiversity monitoring and clarify that the authors’ contribution extends this work through a disease-ecology lens.

#response[We have added reference to these (and additional) monitoring networks in the introduction and highlighted how these can contribute to biosurveillance and provide the basis for deeper integration of biodiversity data in long-term biosurveillance programs.]

=== _Methods/Conceptual Framework_

==== Comment 4:

I like the idea of weighting species differently, and similar approaches have
been attempted in other disease systems (e.g., bird communities in the West Nile virus
system or mammal communities in the Lyme disease system). In many cases, species
weighting decisions are made somewhat ad hoc by modelers. 

I appreciated the ideas introduced in lines 181-188. After reading these lines, I felt that a deeper methodological discussion would be valuable (separate from Box 2 or the supplement example of weight sensitivity analysis). Because this manuscript is presented as a protocol paper, I would welcome a discussion of how species weighting could be implemented in practice, including potential pitfalls and non-negotiable assumptions.

#response[In response to both this comment and a comment from the other reviewer, we have added additional detail on the process of selecting weights, and accounting for uncertainty in weight selection in this part of the manuscript.]

=== _Case Studies_

==== Comment 5:
I recommend ensuring consistency in tone, structure, and figure referencing
across both case studies. At present, they read somewhat as though there were written by
different authors. I particularly appreciated how Case Study 1 clearly walked through how
dose and uncertainty could be interpreted across regions.

#response[We have added more detail focused on map interpretation in order to unify the tone across case studies.]

=== _Discussion_

==== Comment 6: 

The Discussion section feels somewhat rushed compared with the rest of the
manuscript. There are several long stretches without citations, and some examples appear
somewhat selectively chosen. I would also encourage more discussion of funding and
logistical limitations. The Introduction discussion international treaties and policy in some detail, but these ideas are less developed in the Discussion. Aligning the tone of the
Discussion more closely with the Introduction would improve coherence. The final
paragraph in particular could benefit from additional development.

#response[We have expanded upon the discussion and added several citations. This includes better placing our work in the context of the relationship between biogeography and pathogen biology, and further details on the funding model and integration with long-term biodiversity monitoring programs.]

== Minor Comments (Line feedback)

- Line 31: Add a space after punctuation before “We”.

#response[We have implemented this change.]

- Line 37: A citation could be included here. There are people in the wildlife-disease monitoring space that have tried to integrate biodiversity data to guide sampling efforts (e.g., Raina Plowright, Christine Johnson, Cara Brooks, etc).

#response[We have added citations to @Plowright2017PathwaysZoonotic @Johnson2020GlobalShifts and @Eby2023PathogenSpillover after the first clause of this sentence.]

- Line 39: Biodiversity increases can also influence pathogen presence. I would simplify this sentence to refer to biodiversity itself and changes in biodiversity indices.

#response[We have adjusted the phrasing to refer to both the amount of biodiversity and its change.]

- Lines 45-47: Modifying geographic ranges does not always increase cross-species transmission, it can also reduce species ranges or abundances. Consider changing to“potentially leading to increased cross-species transmission”.

#response[We have implemented this change.]

- Lines 59-61: If a distinction is being made between biodiversity surveillance and host ecology, it may help to define these terms earlier (e.g., biodiversity surveillance as species presence records, host ecology as relationships among species and landscapes).

#response[We have refactored this sentence to say "wildlife disease surveillance" instead of "biodiversity surveillance", as that aligns more with the original intended meaning.]

- Lines 63-64: There has been substantial work establishing biodiversity monitoring systems. I encourage the authors to engage more broadly with the literature beyond the current citation network.

#response[We have added additional discussion of long-term biodiversity monitoring programs, e.g. LTER, NEON, BBS, and TERN.]

- Lines 69-71: Related to earlier comments about BON, I suggest revising to “how a long-term biodiversity network perspective would be implemented within the BON system”

#response[We have refactored this sentence to say "how biosurveillance can be implemented within the BON perspective"]

- Lines 77-80: The sentence could be simplified. For example, “the type of data useful for a particularly context depends on existing data and local drivers of pathogen prevalence. For example, historic occurrence records can contextualize contemporary species records, while pathogen prevalence may be informed by the proportion of susceptible hosts”

#response[We have refactored this and the following sentence for clarity.]

- Lines 80-83: The seems to be a central argument of the paper but becomes lost in long sentences. Consider simplifying the logic around trade-offs between occurrence data and pathogen prevalence data.

#response[We have refactored this sentence for clarity.]

- Line 97: Consider replacing “conceptualizing” with “identifying” or “estimating”

#response[We have changed "conceptualizing" to "mitigating" as it aligns closer with our intended meaning.]

- Lines 97-99: “In a state conducive to transmission to humans” feels somewhat redundant. Consider simplifying.

#response[We have removed this clause as transmissability to humans is directly stated in the next clause.]

- Lines 100-101: If discussion local hazard fluctuations, “shifts in host abundance” may be clearer than “host distribution”

#response[We have implemented this change.]

- Lines 106-108: Specify that exposure is not static across space and time

#response[We have stated that these drivers all fluctuate through time.]

- Lines 109-111: To emphasize anthropogenic changes, consider rearranging the sentence to begin with urban encroachment or climate-driven species migration

#response[We have implemented this change.]

- Lines 113-121: While the cited papers are nice and relevant, many originate from the same research group. Including additional perspectives would strengthen the balance of the section.

#response[We have added several new references to this section.]

- Line 123: Some more recent literature could be cited here. Additionally, because community composition does not always influence pathogen dynamics, the statement may need to be softened or balanced with additional citations.

#response[We have added citations from #cite(<Johnson2013BiodiversityDecreases>, form: "prose")#cite(<Civitello2015BiodiversityInhibits>, form: "prose"), #cite(<Gibb2020ZoonoticHost>, form: "prose") and #cite(<Keesing2021ImpactsBiodiversity>, form: "prose") and adjusted the language to say "host communities _can be_ a central determinant of pathogen dynamics..."]

- Lines 128-129: I particularly enjoyed this sentence

#response[We thank the reviewer for their comment.]

- Lines 130-131: Because this introduces the methods, consider breaking it into shorter sentences for clarity

#response[We have implemented this change.]


- Lines 137-138: Clarify that when pathogen prevalence relationships are known but abundance is not, sampling should prioritize host identity.

#response[We are a bit confused about what the reviewer means here. If pathogen prevalence is known, wouldn't host identity already be known because prevalence data reflects the proportion of individuals of a given host species in which the pathogen was found? Does the reviewer mean host identity should be prioritized for other hosts outside of those with available prevalence data? Or that the abundance of the host with available prevalence data should be the next priority for sampling?

We have not implemented any changes here, but would be willing to if the reviewer what clarify their comment.]

- Lines 144-145: “Quantitative” may be redundant, consider removing or replacing with
“statistical”

#response[We have implemented this change.]

- Line 159: Adding a concrete example (e.g., overlapping susceptible host populations) may help

#response[We have added an example of red deer as a "bridge host" of bovine tuberculosis, citing @Nugent2011MaintenanceSpillover, to emphasize this point.]

- Lines 161-163: Pathogen extent & host susceptibility weighting seem conceptually distinct and could be clarified

#response[We have restructured this sentence for clarity.]

- Line 168: simplify to “both across space and in the type of data collected”
- Line 169: Replace “this” with “our protocol”

#response[We have implemented this change.]

- Lines 184-185: I liked the idea of weighting host competency based on testing protocol (eg serology vs PCR)

#response[We thank the reviewer for their comment, and have added additional details here in response to a comment from the other reviewer.]

- Line 200: refer simply to figure 1 if panel order is not sequential

#response[We have implemented this change.]

- Line 201-204: These important points might be introduced earlier in the Introduction

#response[We have added a sentence in the introduction discussing how the protocol enables prioritization of the type of data to collect at each site.]

- Lines 216-221: it may be useful to mention that pathogen prevalence and host distribution are often reported at different spatial scales. For example pathogen prevelance may be reported as a lat/lon, but that may represent the centroid of a plot (with unknown size) or could be at the regional administrative level.

#response[We now mention this at the beginning of the subsequent paragraph.]

- Line 258: Consider removing “(which is maintained by the first author)” or moving that detail to the Data Availability section

#response[We have implemented this change.]

- Lines 300-303: As a preference, it might be helpful to present the more data-rich case study (south korea) first

#response[We have considered the reviewers suggestion, but we have decided against this because leading with the less data rich case study first enables introduction of dose/uncertainty maps first, without the additional complexity contributed by accounting for prevalence uncertainty in the second case study. An attempt to restructure them left the India case study feeling a bit anticlimatic in our opinion, so we have kept the structure of India first and Korea second. ]

- Lines 308-309: if a package or specific download citation was used to retrieve GBIF data, it
should be cited (often GBIF gives you a specific citation based on the occurrence data
downloaded)

#response[We have included the GBIF DOIs associated with each study here.]

- Lines 326-331: These sentences could likely be written more concisely

#response[We have refactored these sentences for clarity.]

- Line 340: replace “this” with “the map of India”

#response[We have implemented this change.]

- Lines 340-343: the map is well explained, but a simpler bivariate legend (e.g., 3x3) might help readers quickly interpret patterns. Supplemental figures could provide alternative visualizations (I get that it’s nice to have 50% cutoff for sampling location priority)

#response[We have adjusted the dose/uncertainty bivariate maps for each case study to use a 3x3 format.]

- Lines 343-345: if “predicted suitability” corresponds to dose, consider stating this explicitly

#response[We have explicitly stated this is biodiversity dose for clarity.]

- Lines 347-350: Overlaying the sampling points in teal was very helpful.

#response[We thank the reviewer for their comment.]

- Lines 363-366: it may help to clarify why these particular papers were chosen. Did they only look at your specific host-virus pairs or did they represent certain geographic regions?

#response[We have added that these studies are selected because they are the sources of data in the ArHa database in South Korea with explict prevalence estimates.]

- Lines 407-416: I appreciated that temporal limitations were addressed early in the Discussion

#response[We thank the reviewer for their comment.]

- Lines 429-430: As noted earlier, many long-term monitoring programs already exist. The manuscript could frame these are opportunities to extend biodiversity monitoring toward biosurveillance in a more standardized format.

#response[We have adjusted this sentence to emphasize the protocol enables kickstarting monitoring in the low data context, and integration with existing biodiversity monitoring programs.]


- Lines 438-440: I encourage the authors to be cautious with claims that this work represents the “first step” towards biosurveillance. Similar ideas and regional analyses already exist in the literature, even if implemented with different tools or datasets.

#response[We agree with the reviewer and have rephrased this as a necessary next step for integrating biosurveillance with biodiversity monitoring.]

- Line 442: Minor figure design suggestion, aligning the arrows between panels (the tips of arrows directly touch) may improve visual clarity

#response[We have adjusted the arrow pointing from species weighting to the path pointing to dose and dose uncertainty for clarity.]

- Line 448: Spell out “and” instead of using “+” . Also include citations for country boundary datasets if applicable.

#response[We have implemented this change.]

- Line 459: Box 1 for SDMs is nicely written and will be easy for people to point trainees in that
direction for a gentle overview on the method.

#response[We thank the reviewer for their comment.]

- Line 568: Include citations for country boundary datasets if applicable.

#response[We have included that GADM is the source of boundaries in both case study figure captions.]


= Reviewer 2

Developing surveillance systems for wildlife pathogens continues to be a tricky problem in the field of One Health due to the presence of many unknowns and the unavoidable constraints in sampling. One cannot sample all species all at once, so there must be some prioritization system. The authors of this manuscript have proposed a method for designing wildlife pathogen surveillance systems that stakeholders could use to better target some of the uncertainties about hosts or pathogen prevalence in a given area. The protocol appears to be ideal for use in an evolving information environment, where new host or prevalence can be added to the underlying model to create the next iteration of sampling prioritization. The manuscript followed a somewhat unusual format for a research article, but it worked well for how the content was presented, as an introduction to the concepts and the model, then a demonstration using two case studies. The boxes provided some additional useful information about how to interpret species distribution models and the recommended sampling locations outputted from the model. Overall, I think this protocol could be an excellent addition to the field. The manuscript text was clear and the figures illustrated some of the key outputs and concepts discussed in the text, such that users could refer back to this manuscript to help with interpretation of their own output. I did not have any major concerns about the model, though I did think there were many places where some additional explanation about how choices are made about species weights and other model variables would be valuable for potential users.

#response[We thank the reviewer for the kind assessment.]

One other general comment: would it be useful to have a glossary of terms somewhere in the manuscript or supplement? There are some esoteric terms that get used throughout the paper and are defined in the text, but a centralized location for these terms could help readers with comprehension. Key terms may include: risk, hazard, spillover, biosurveillance, discovery sampling, prevalence sampling, Biodiversity Observation Network, biodiversity dose.

#response[We have added a glossary in the supplemental material of the manuscript.]

L19-20: Expand on this statement a little to provide more justification. As it stands now, this seems like a fairly empty proclamation. However, the first paragraph of Section 2 does a fantastic job with distinguishing hazard from risk and justifying why biosurveillance is a key part of estimating hazard. Pulling in some of that nuance here would make a stronger case.

#response[We have added that this specifically is aimed to identify regions whwere spillover is most probable to target preventative measures.]

L32: Hyphen not needed in Arenaviridae. Also make sure to italicize virus family names. Lastly, and this may be pedantic, but “small mammals” or “rodents and shrews” would be better than rodents here and throughout the paper. Suncus murinus and Crocidura lasiura are shrews (Eulipotyphla), not rodents. And I think broadly the hantaviruses can infect both taxa, so you should be as inclusive as possible, even if arenaviruses are more specific to rodents.

#response[We've adjusted both family names and changed to using both "small mammals" and "rodents and shrews" throughout the manuscript, rather than strictly "rodents".]

L34-35: This is more of a general comment, but it was not immediately clear how this protocol would get used in the public health space. Who are the intended users of this protocols and how would the information about sampling and testing be distributed to stakeholders? Do you think this is a purely academic exercise, or is it intended that governments and agencies would be leading these efforts? It may be worth some additional discussion at the end of the paper to make the strong case to specific end users.

#response[We have added an additional paragraph in the discussion about pragmatic details of funding and integration with existing biodiversity monitoring programs, in conjunction with comments from the other reviewer.]

L124-127: Why is this considered part of risk? As laid out in the sections above, this seems more in line with hazard rather than risk, because risk is largely driven by anthropogenic activities, no? Changes in the host side of things seem like they would have more influence on hazard than risk.

#response[We have adjusted this sentence to explicitly say the density of infected hosts in a location proximate to human populations, to make it explicit that why we are talking about risk rather than hazard.]

L184-185: Unclear which of the methods is the one with higher false-negatives. One would initially guess that PCR has the higher false-negative rate because of the frequently short window of pathogen shedding, such that it may be difficult to find a host actively shedding, even if they are ultimately a competent host. Whereas with serology, one runs the risk of including many noncompetent hosts that have simply been exposed but the pathogen does not replicate in or transmit from this host. So those would be false-positives in the eyes of many ecologists. This question seems like it would also have a fairly large influence on the model uncertainty. Is this something that needs to be included in sensitivity analyses? And is there a way of incorporating uncertainty in host weights into the calculation of biodiversity dose, beyond just random noise, but a more specific way of quantifying differences in uncertainty among host species?

#response[We have removed the reference to one method having higher false-negatives and added several sentences of additional detail to clarify the difference in what serological vs. PCA data are saying, and more on what the implications of each type of data for sampling prioritization are.]


L199: What about prevalence uncertainty? If a small sample at a location showed high prevalence, it does not mean it is always that high and is likely biased. This location would have high prevalence uncertainty (10/1000 has less uncertainty than 1/100). Is there any way to include prevalence uncertainty into this model, e.g., weighting inversely by sample size? If it is already included in the model, a short section explaining how it is done would be helpful.

#response[We agree with the reviewer that this is a highly relevant point, which we demonstrate in the second case study. We have added an additional sentence on this in the "Integration of prevalence data into monitoring network design", and it is mentioned explicitly in Supplemental Table 1.]

L213: Could the authors provide more guidance about how all these uncertainties about host competence and prevalence get boiled down into w\_{i}? Is there an equation that could be built to estimate w\_{i}, or are the values based solely on investigator judgement?

#response[We have added additional discussion in the "Biodiversity dose, host weighting, and uncertainty" section of the protocol to address this question, specifically focusing on ways to assess weight sensitivity to larger, biologically meaningful differences in weight selection and host inclusion/exclusion. We have added a recommendation for future work on assessing explicit formulae for constructing weights based on various forms of available data (e.g. phylogenies, serology data, and the capacity of species to serve as hosts) using simulation as a sandbox to test the efficacy of different formulae for weight construction, although we believe this is out of scope of this manuscript.]

L281-283: Could you also apply this to investigate the sensitivity of maps to changes in weights that do indeed reflect biological differences in host relevance, such as when one doesn’t know whether all species are equally competent hosts? Would there need to be something larger than small “nudges” to the weights?

#response[ 
We agree with the reviewer that other forms of sensitivity analysis are prudent to assess the inclusion/exclusion of particular hosts, and larger differences in selected weights that reflect meaningful biological difference in host capacity.

The "nudge" version of sensitivity analysis is mostly to ensure the selected weights are not near a "steep" region in weight-space where the output changes dramatically based.

We have added additional detail in the "Biodiversity dose, host weighting, and uncertainty" section of the manuscript to address how uncertainty in weight selection can be account for, namely by creating priority maps for each set of competing weights and included hosts, and targeting sampling toward the regions that are high priority under each set of candidate priority maps.  
]

L300-L303: Maybe indicate that the South Korea case study is about hosts of Hantaviridae. It’s not currently stated explicitly here.

#response[We have implemented this change.]

L308-309: If there is no prevalence or sampling data in some locations, how can one be sure about which hosts are contributing to pathogen maintenance? Are there potentially unknown host species that need to be accounted for? And how would one go about assigning weights to these species?

#response[We have added several sentences about the choices made in host inclusion and uncertainty about weight selection in the "Biodiversity dose, host weighting, and uncertainty" section of the protocol that we believe addresses these concerns with more methodological detail. ]

L335-338: I understand what is trying to be done here with this case study, but it sort of falls short of being a great example of how one might approach this type of low information sampling. This seems from a glance to be just an assortment of relatively common species in India. The authors might have to provide a disclaimer that the species chosen by Rodrigues may be different than what would have been chosen by an alternative method that uses all available biodiversity information and phylogenetic proximity to known hosts.

#response[We agree that with modern available data there would be several additional potential hosts we would include, and have added a disclaimer to emphasize we only use these species for the sake of example because they are the species sampled by Rodrigues et al. and a contemporary approach would use addition information to select hosts.]

L338: What were the chosen weights? Which species ended up with the highest weight in the images shown in Figure 2?

#response[We have added the values of the randomly drawn weights in the text.]

L348-349: What is the process for deciding the tilting value alpha? What sort of behavior in the model is “ideal” and how would a user decide?

#response[We have added a sentence on how this corresponds to mean/variance of priority values and related this to uncertainty in decisions about what hosts to include and their selected weights, and included an additional figure in the Supplemental Material that visualizes the distribution of priority value across selected sites for different values of the tilting parameter $alpha$.]

L359: What is the justification for including Imjin virus here? Is there any evidence that this virus can or does infect humans? Or is this operating from a pure hazard viewpoint that any hantavirus in a host can potentially cause human infections? And what do we know about the role of these species in transmission? Are they the only hosts, or just the dominant hosts? How many additional species contribute to hantavirus transmission in South Korea that are not included in this case study?

#response[We have included Imjin virus because there was good existing prevalence data, and there is evidence that it can replicate in human macrophages @Shin2012DistinctInnate. The goal was to show a sampling prioritization that reflects a balance of sampling known zoonoses and viruses of potential spillover concern, and we have edited the text to include this justification.]

L363-366: If the authors were performing this case study as a real analysis or sample planning exercise, what sort of sensitivity analyses would the authors suggest running to incorporate additional known hosts and other potential hosts with unknown status?

#response[We have added additional discussion regarding other forms of sensitivity analysis to assess the impact of host inclusion/exclusion in the "Biodiversity dose, host weighting, and uncertainty" section of the protocol, namely comparing priority maps with different species included and excluded and targeting sampling toward the regions that are high priority under all host pools.]

L384: Similar comment as above, how was the tilting value alpha chosen for this case study?

#response[We have added a parathetical to refer to the newly added comment on tilting value selection in the India case study here.]

L481-483: The way I thought about this was in terms of additive probabilities. Because the presence of species are not mutually exclusive, their additive probabilities would be P(A or B) = P(A) + P(B) – P(A and B). However, since P(A and B) requires a joint SDM to estimate, it does not get subtracted off in the probability estimate and therefore we get an overestimate of sympatry between two species. And obviously this gets compounded when you stack an even larger number of potential sympatric species.

#response[Interestingly, overall species richness is often better predicted by Stacked SDMs instead of Joint SDMs, and the difference between methods is typically smaller than the difference induced by using different stacking procedures (e.g. thresholding to create binary ranges vs. using continuous suitability scores) as found by #cite(<Zurell2020TestingSpecies>).

We have added this caveat and reference to this section of the manuscript.
]


#bibliography("refs.bib", style: "springer-basic-author-date")