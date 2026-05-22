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

#let add(body) = text(fill: rgb("#1b8e1b"), weight: "regular")[#body]
#let change(body) = text(fill: rgb("#2b2ba5"))[#underline(body, stroke: rgb(0, 90, 90))]
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

#response[We thank the reviewer for their kind words summarizing the manuscript's strengths.]

I was able to clone the repository and run part of the code. Without troubleshooting or
actively fixing errors, I was able to run the India case study but not the South Korea case
study. I initially suspected the issue might relate to OSM object, but it appears to be
unrelated. However, I am primarily an R user, so this could reflect my own skillset rather
than an issue with the repository itself. If may be helpful for the authors to ask another
primarily R user to test the workflow on a different machine.

#response[*TODO* BAR. We believe that this is related to an update that occurred from Julia 1.11 (which is where the original code was developed) and Julia 1.12. ]

== Major Comments (General Feedback)

Overall, I found this article to be well written, visually appealing, and clearly motivated by policy relevance in the Introduction. 

However, I repeatedly found myself getting stuck on the framing that BON – and this protocol in particularly – is the primary or only framework addressing these questions. There has been substantial effort & investment in standardized biodiversity sampling, and there have also been many attempts to link biodiversity data with disease ecology indices (e.g., pathogen prevalence, pathogen load, host community structure). 

Because the manuscript does not sufficiently situate itself within the broader literature or acknowledge earlier related work (outside of the co-author list), the framing occasionally felt somewhat narrow. This was frustrating at time because in most other respects the manuscript is quite strong. 

I recognize that several authors on this paper have made important contribution to the biodiversity-disease space (and I enjoy these papers!). However, this is a large and active field with many research groups working on related problems. I would encourage the authors to broaden the citation base so that the manuscript better reflects that broader community. A brief scan of the reference list suggests a substantial number of citations include members of the author team (even with the “et al” truncation). 

#response[TODO we have added additional citations from]

Similarly, the manuscript highlights BON as the primary monitoring framework, despite the existence of many other long-term biodiversity monitoring networks. While BON is an important example, it would be more appropriate to frame it as one of several efforts working toward standardized biodiversity monitoring. This work builds on substantial existing foundation, and acknowledging the broader literature would strengthen the manuscript.

== Specific Comments

=== _Introduction_

==== Comment 1:

More consideration could be given to citation choices in the opening paragraph of the Introduction. While the cited research groups have contributed significantly to this field, it would be beneficial to diversify the research groups cited in the opening paragraph, given the breadth of work in biodiversity-disease research.

==== Comment 2: 

After reading lines 49-57, I wondered who would realistically fund the surveillance and standardization systems being proposed. Are these cost estimates or funding mechanisms that could support implementation of the monitoring protocol suggested here? While these logistical considerations may be too detailed for the Introduction, I hope the Discussion includes a more substantial consideration of the trade- offs between ideal surveillance systems and what is practical to implement. Treaties and global initiatives may call for these systems, but the protocol proposed here implies substantial additional effort (at the individual and institutional levels), and it would strengthen the paper to acknowledge the logistical challenges users may face.

#response[TODO: we add a paragraph in the discussion on this and relate to financing GBiOS and GBF stuff]

==== Comment 3: 

The manuscript focuses heavily on BON, but there are many long-term
monitoring networks conducting similar work. It would be useful to frame this approach as
building on existing monitoring infrastructure rather than positioning BON as the primary
monitoring network. For example, within North America alone there are several well-
established networks including US LTER, US NEON, North American Breeding Bird Survey.
Highlighting similar networks would help acknowledge the substantial investments that
have already been made in biodiversity monitoring and clarify that the authors’ contribution extends this work through a disease-ecology lens.

#response[TODO: we acknowledge much work has been done outside the BON perspective and add this to the intro and relevant places]

=== _Methods/Conceptual Framework_

==== Comment 4:

I like the idea of weighting species differently, and similar approaches have
been attempted in other disease systems (e.g., bird communities in the West Nile virus
system or mammal communities in the Lyme disease system). In many cases, species
weighting decisions are made somewhat ad hoc by modelers. 

#response[TODO: we thank the reviewer for pointing us to this work and we have added reference to this in ...]

I appreciated the ideas introduced in lines 181-188. After reading these lines, I felt that a deeper methodological discussion would be valuable (separate from Box 2 or the supplement example of weight sensitivity analysis). Because this manuscript is presented as a protocol paper, I would welcome a discussion of how species weighting could be implemented in practice, including potential pitfalls and non-negotiable assumptions.

#response[TODO: yeah we should do this]

=== _Case Studies_

==== Comment 5:
I recommend ensuring consistency in tone, structure, and figure referencing
across both case studies. At present, they read somewhat as though there were written by
different authors. I particularly appreciated how Case Study 1 clearly walked through how
dose and uncertainty could be interpreted across regions.

#response[TODO: we have made them read more like the same person.]

=== _Discussion_

==== Comment 6: 

The Discussion section feels somewhat rushed compared with the rest of the
manuscript. There are several long stretches without citations, and some examples appear
somewhat selectively chosen. I would also encourage more discussion of funding and
logistical limitations. The Introduction discussion international treaties and policy in some detail, but these ideas are less developed in the Discussion. Aligning the tone of the
Discussion more closely with the Introduction would improve coherence. The final
paragraph in particular could benefit from additional development.

#response[TODO: we discuss funding and adjusted the tone to match the introduction more closely]

== Minor Comments (Line feedback)

- Line 31: Add a space after punctuation before “We”.

- Line 37: A citation could be included here. There are people in the wildlife-disease monitoring space that have tried to integrate biodiversity data to guide sampling efforts (e.g., Raina Plowright, Christine Johnson, Cara Brooks, etc).

- Line 39: Biodiversity increases can also influence pathogen presence. I would simplify this sentence to refer to biodiversity itself and changes in biodiversity indices.

- Lines 45-47: Modifying geographic ranges does not always increase cross-species transmission, it can also reduce species ranges or abundances. Consider changing to“potentially leading to increased cross-species transmission”.

- Lines 59-61: If a distinction is being made between biodiversity surveillance and host ecology, it may help to define these terms earlier (e.g., biodiversity surveillance as species presence records, host ecology as relationships among species and landscapes).

- Lines 63-64: There has been substantial work establishing biodiversity monitoring systems. I encourage the authors to engage more broadly with the literature beyond the current citation network.

- Lines 69-71: Related to earlier comments about BON, I suggest revising to “how a long-term biodiversity network perspective would be implemented within the BON system”

- Lines 77-80: The sentence could be simplified. For example, “the type of data useful for a particularly context depends on existing data and local drivers of pathogen prevalence. For example, historic occurrence records can contextualize contemporary species records, while pathogen prevalence may be informed by the proportion of susceptible hosts”

- Lines 80-83: The seems to be a central argument of the paper but becomes lost in long sentences. Consider simplifying the logic around trade-offs between occurrence data and pathogen prevalence data.

- Line 97: Consider replacing “conceptualizing” with “identifying” or “estimating”

- Lines 97-99: “In a state conducive to transmission to humans” feels somewhat redundant. Consider simplifying.

- Lines 100-101: If discussion local hazard fluctuations, “shifts in host abundance” may be clearer than “host distribution”

- Lines 106-108: Specify that exposure is not static across space and time

- Lines 109-111: To emphasize anthropogenic changes, consider rearranging the sentence to begin with urban encroachment or climate-driven species migration

- Lines 113-121: While the cited papers are nice and relevant, many originate from the same research group. Including additional perspectives would strengthen the balance of the section.

- Line 123: Some more recent literature could be cited here. Additionally, because community composition does not always influence pathogen dynamics, the statement may need to be softened or balanced with additional citations.

- Lines 128-129: I particularly enjoyed this sentence

- Lines 130-131: Because this introduces the methods, consider breaking it into shorter sentences for clarity

- Lines 137-138: Clarify that when pathogen prevalence relationships are known but abundance is not, sampling should prioritize host identity.

- Lines 144-145: “Quantitative” may be redundant, consider removing or replacing with
“statistical”

- Line 159: Adding a concrete example (e.g., overlapping susceptible host populations) may help

- Lines 161-163: Pathogen extent & host susceptibility weighting seem conceptually distinct and could be clarified

- Line 168: simplify to “both across space and in the type of data collected”
- Line 169: Replace “this” with “our protocol”
- Lines 184-185: I liked the idea of weighting host competency based on testing protocol (eg serology vs PCR)

- Line 200: refer simply to figure 1 if panel order is not sequential

- Line 201-204: These important points might be introduced earlier in the Introduction

- Lines 216-221: it may be useful to mention that pathogen prevalence and host distribution are often reported at different spatial scales. For example pathogen prevelance may be reported as a lat/lon, but that may represent the centroid of a plot (with unknown size) or could be at the regional administrative level.

- Line 258: Consider removing “(which is maintained by the first author)” or moving that detail to the Data Availability section

- Lines 300-303: As a preference, it might be helpful to present the more data-rich case study (south korea) first

- Lines 308-309: if a package or specific download citation was used to retrieve GBIF data, it
should be cited (often GBIF gives you a specific citation based on the occurrence data
downloaded)

- Lines 326-331: These sentences could likely be written more concisely

- Line 340: replace “this” with “the map of India”

- Lines 340-343: the map is well explained, but a simpler bivariate legend (e.g., 3x3) might help readers quickly interpret patterns. Supplemental figures could provide alternative visualizations (I get that it’s nice to have 50% cutoff for sampling location priority)

- Lines 343-345: if “predicted suitability” corresponds to dose, consider stating this explicitly

- Lines 347-350: Overlaying the sampling points in teal was very helpful.

- Lines 363-366: it may help to clarify why these particular papers were chosen. Did they only look at your specific host-virus pairs or did they represent certain geographic regions?

- Lines 407-416: I appreciated that temporal limitations were addressed early in the Discussion

- Lines 429-430: As noted earlier, many long-term monitoring programs already exist. The manuscript could frame these are opportunities to extend biodiversity monitoring toward biosurveillance in a more standardized format.

- Lines 438-440: I encourage the authors to be cautious with claims that this work represents the “first step” towards biosurveillance. Similar ideas and regional analyses already exist in the literature, even if implemented with different tools or datasets.

- Line 442: Minor figure design suggestion, aligning the arrows between panels (the tips of arrows directly touch) may improve visual clarity

- Line 448: Spell out “and” instead of using “+” . Also include citations for country boundary datasets if applicable.

- Line 459: Box 1 for SDMs is nicely written and will be easy for people to point trainees in that
direction for a gentle overview on the method.

- Line 568: Include citations for country boundary datasets if applicable.

= Reviewer 2

Developing surveillance systems for wildlife pathogens continues to be a tricky problem in the field of One Health due to the presence of many unknowns and the unavoidable constraints in sampling. One cannot sample all species all at once, so there must be some prioritization system. The authors of this manuscript have proposed a method for designing wildlife pathogen surveillance systems that stakeholders could use to better target some of the uncertainties about hosts or pathogen prevalence in a given area. The protocol appears to be ideal for use in an evolving information environment, where new host or prevalence can be added to the underlying model to create the next iteration of sampling prioritization. The manuscript followed a somewhat unusual format for a research article, but it worked well for how the content was presented, as an introduction to the concepts and the model, then a demonstration using two case studies. The boxes provided some additional useful information about how to interpret species distribution models and the recommended sampling locations outputted from the model. Overall, I think this protocol could be an excellent addition to the field. The manuscript text was clear and the figures illustrated some of the key outputs and concepts discussed in the text, such that users could refer back to this manuscript to help with interpretation of their own output. I did not have any major concerns about the model, though I did think there were many places where some additional explanation about how choices are made about species weights and other model variables would be valuable for potential users.

#response[We thank the reviewer for the kind assessment.]

One other general comment: would it be useful to have a glossary of terms somewhere in the manuscript or supplement? There are some esoteric terms that get used throughout the paper and are defined in the text, but a centralized location for these terms could help readers with comprehension. Key terms may include: risk, hazard, spillover, biosurveillance, discovery sampling, prevalence sampling, Biodiversity Observation Network, biodiversity dose.

#response[TODO: we might do this idk]

L19-20: Expand on this statement a little to provide more justification. As it stands now, this seems like a fairly empty proclamation. However, the first paragraph of Section 2 does a fantastic job with distinguishing hazard from risk and justifying why biosurveillance is a key part of estimating hazard. Pulling in some of that nuance here would make a stronger case.

L32: Hyphen not needed in Arenaviridae. Also make sure to italicize virus family names. Lastly, and this may be pedantic, but “small mammals” or “rodents and shrews” would be better than rodents here and throughout the paper. Suncus murinus and Crocidura lasiura are shrews (Eulipotyphla), not rodents. And I think broadly the hantaviruses can infect both taxa, so you should be as inclusive as possible, even if arenaviruses are more specific to rodents.

#response[TK: we've adjusted both family names and changed to using both "small mammals" and "rodents and shrews" throughout the manuscript, rather than strictly "rodents"]

L34-35: This is more of a general comment, but it was not immediately clear how this protocol would get used in the public health space. Who are the intended users of this protocols and how would the information about sampling and testing be distributed to stakeholders? Do you think this is a purely academic exercise, or is it intended that governments and agencies would be leading these efforts? It may be worth some additional discussion at the end of the paper to make the strong case to specific end users.

#response[TK: we have added discussion points about pragmatic constraints in lines X - X, in conjunction with comments from the other reviewer]

L124-127: Why is this considered part of risk? As laid out in the sections above, this seems more in line with hazard rather than risk, because risk is largely driven by anthropogenic activities, no? Changes in the host side of things seem like they would have more influence on hazard than risk.

L184-185: Unclear which of the methods is the one with higher false-negatives. One would initially guess that PCR has the higher false-negative rate because of the frequently short window of pathogen shedding, such that it may be difficult to find a host actively shedding, even if they are ultimately a competent host. Whereas with serology, one runs the risk of including many noncompetent hosts that have simply been exposed but the pathogen does not replicate in or transmit from this host. So those would be false-positives in the eyes of many ecologists. This question seems like it would also have a fairly large influence on the model uncertainty. Is this something that needs to be included in sensitivity analyses? And is there a way of incorporating uncertainty in host weights into the calculation of biodiversity dose, beyond just random noise, but a more specific way of quantifying differences in uncertainty among host species?

#response[TK: long story short.]


L199: What about prevalence uncertainty? If a small sample at a location showed high prevalence, it does not mean it is always that high and is likely biased. This location would have high prevalence uncertainty (10/1000 has less uncertainty than 1/100). Is there any way to include prevalence uncertainty into this model, e.g., weighting inversely by sample size? If it is already included in the model, a short section explaining how it is done would be helpful.

#response[TK This is a good point that is mentioned in the supplemental table. We have added an explicit mention of the uncertainty in prevalence estimates due to sample size in lines TK of the main text.]

L213: Could the authors provide more guidance about how all these uncertainties about host competence and prevalence get boiled down into w_{i}? Is there an equation that could be built to estimate w_{i}, or are the values based solely on investigator judgement?

#response[TK. ]

L281-283: Could you also apply this to investigate the sensitivity of maps to changes in weights that do indeed reflect biological differences in host relevance, such as when one doesn’t know whether all species are equally competent hosts? Would there need to be something larger than small “nudges” to the weights?

#response[TK. The nudge version of sensitivity analysis is mostly to ensure the selected weights are not near a ]


L300-L303: Maybe indicate that the South Korea case study is about hosts of Hantaviridae. It’s not currently stated explicitly here.

#response[TK we did this]

L308-309: If there is no prevalence or sampling data in some locations, how can one be sure about which hosts are contributing to pathogen maintenance? Are there potentially unknown host species that need to be accounted for? And how would one go about assigning weights to these species?

L335-338: I understand what is trying to be done here with this case study, but it sort of falls short of being a great example of how one might approach this type of low information sampling. This seems from a glance to be just an assortment of relatively common species in India. The authors might have to provide a disclaimer that the species chosen by Rodrigues may be different than what would have been chosen by an alternative method that uses all available biodiversity information and phylogenetic proximity to known hosts.

#response[TODO: we agree and have acknowledged contemporary approach would consider different species]

L338: What were the chosen weights? Which species ended up with the highest weight in the images shown in Figure 2?

#response[TODO: we have added the values of the randomly drawn weights in the text]

L348-349: What is the process for deciding the tilting value alpha? What sort of behavior in the model is “ideal” and how would a user decide?

#response[TODO: add a paragraph on what tilting really means: look at the distribution of values in feature and uncertainty space]

L359: What is the justification for including Imjin virus here? Is there any evidence that this virus can or does infect humans? Or is this operating from a pure hazard viewpoint that any hantavirus in a host can potentially cause human infections? And what do we know about the role of these species in transmission? Are they the only hosts, or just the dominant hosts? How many additional species contribute to hantavirus transmission in South Korea that are not included in this case study?

L363-366: If the authors were performing this case study as a real analysis or sample planning exercise, what sort of sensitivity analyses would the authors suggest running to incorporate additional known hosts and other potential hosts with unknown status?

L384: Similar comment as above, how was the tilting value alpha chosen for this case study?

L481-483: The way I thought about this was in terms of additive probabilities. Because the presence of species are not mutually exclusive, their additive probabilities would be P(A or B) = P(A) + P(B) – P(A and B). However, since P(A and B) requires a joint SDM to estimate, it does not get subtracted off in the probability estimate and therefore we get an overestimate of sympatry between two species. And obviously this gets compounded when you stack an even larger number of potential sympatric species.

#response[This is why its importance to remember that this is a suitability score and not a probability of occurrence. TK: add back cite of Zurell paper that shows SSDMs are better at richness estimates]

= Editor



#enum()[
  *Please amend your detailed Financial Disclosure statement. This is published with the article. It must therefore be completed in full sentences and contain the exact wording you wish to be published.*
  #enum(numbering: "i.")[
    Please clarify all sources of financial support for your study. List the grants, grant numbers, and organizations that funded your study, including funding received from your institution. Please note that suppliers of material support, including research materials, should be recognized in the Acknowledgements section rather than in the Financial Disclosure. 
  ][
    State the initials, alongside each funding source, of each author to receive each grant. For example: "This work was supported by the National Institutes of Health (\#\#\# to AM; \#\#\# to CJ) and the National Science Foundation "\#\#\#" to AM."
  ][State what role the funders took in the study. If the funders had no role in your study, please state: “The funders had no role in study design, data collection and analysis, decision to publish, or preparation of the manuscript.”
  ][
    If any authors received a salary from any of your funders, please state which authors and which funders.
  ]
][
  *We ask that a manuscript source file is provided at Revision. Please upload your manuscript file as a .doc, .docx, .rtf or .tex.*
][
  *Please upload separate figure files in .tif or .eps format. Also, remove the figures from your manuscript file but keep the legends.*
][
  *Please provide an Author Summary. This should appear in your manuscript between the Abstract (if applicable) and the Introduction, and should be 150–200 words long. The aim should be to make your findings accessible to a wide audience that includes both scientists and non-scientists. Sample summaries can be found on our website under Submission Guidelines*
][
  *We have noticed that you have uploaded Supporting Information files, but you have not included a list of legends. Please add a full list of legends for your Supporting Information files after the references list.*
][
  *Some material included in your submission may be copyrighted. According to PLOS’s copyright policy, authors who use figures or other material (e.g., graphics, clipart, maps) from another author or copyright holder must demonstrate or obtain permission to publish this material under the Creative Commons Attribution 4.0 International (CC BY 4.0) License used by PLOS journals. Please closely review the details of PLOS’s copyright requirements here: PLOS Licenses and Copyright. If you need to request permissions from a copyright holder, you may use PLOS's Copyright Content Permission form. \ \ Please respond directly to this email or email the journal office and provide any known details concerning your material's license terms and permissions required for reuse, even if you have not yet obtained copyright permissions or are unsure of your material's copyright compatibility.*

  #enum(numbering: "a.")[
    Figure 1: Please confirm whether you drew the images / clip-art within the figure panels by hand. If you did not draw the images, please provide (a) a link to the source of the images or icons and their license / terms of use; or (b) written permission from the copyright holder to publish the images or icons under our CC-BY 4.0 license. Alternatively, you may replace the images with open source alternatives. See these open source resources you may use to replace images / clip-art:
  ][
    Figures 2, 3, 4 and Figure 1 in supplement.pdf: please (a) provide a direct link to the base layer of the map (i.e., the country or region border shape) and ensure this is also included in the figure legend; and (b) provide a link to the terms of use / license information for the base layer image or shapefile. We cannot publish proprietary or copyrighted maps (e.g. Google Maps, Mapquest) and the terms of use for your map base layer must be compatible with our CC-BY 4.0 license. 
  ][
    Note: if you created the map in a software program like R or ArcGIS, please locate and indicate the source of the basemap shapefile onto which data has been plotted.

    If your map was obtained from a copyrighted source please amend the figure so that the base map used is from an openly available source. Alternatively, please provide explicit written permission from the copyright holder granting you the right to publish the material under our CC-BY 4.0 license.

    Please note that the following CC BY licenses are compatible with PLOS license: CC BY 4.0, CC BY 2.0 and CC BY 3.0, meanwhile such licenses as CC BY-ND 3.0 and others are not compatible due to additional restrictions. 

    If you are unsure whether you can use a map or not, please do reach out and we will be able to help you. The following websites are good examples of where you can source open access or public domain maps: 

    - U.S. Geological Survey (USGS) - All maps are in the public domain. (http://www.usgs.gov) 
    - PlaniGlobe - All maps are published under a Creative Commons license so please cite “PlaniGlobe, http://www.planiglobe.com, CC BY 2.0” in the image credit after the caption. (http://www.planiglobe.com/?lang=enl) 
    - Natural Earth - All maps are public domain. (http://www.naturalearthdata.com/about/terms-of-use/)
  ]

  #response[
    The maps used in Figures 2, 3, and 4, and the maps of Korea for Supp. Fig 1, all use GADM, which explicitly states that publishing them under CC-BY is allowed #text(fill: blue)[#link("https://gadm.org/license.html", "on their site")].

    We have replaced the zoomed-in maps of selected sites in Figure 4 with maps from OpenStreetMap, which are 

    TK : I am adding the zoom in land cover maps via sentinel imagery
  ]
]
