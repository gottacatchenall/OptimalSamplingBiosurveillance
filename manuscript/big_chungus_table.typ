#show table.cell: set text(size: 0.92em)
#show table.cell.where(y: 0): set text(size: 1.1em)
#show table.cell.where(x: 0): set text(size: 1.1em)

#rotate(
  -90deg,
  reflow: true,

  figure(numbering: "1", caption: [Methods for incorporating prevalence data with varying levels of spatial and temporal coverage])[
  #table(
  columns: (8%, auto, auto, auto, auto),
  rows: (7%, auto, auto, auto, auto),
  inset: 12pt,
  align: (x, y) => {
    if x > 0 and y > 0 {
      left 
    }
    else {
      horizon 
    }
  },
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
  table.header(
    [], 
    [Low Predicted Prevalence \ + Low Prevalence Uncertainty], 
    [Low Predicted Prevalence \ + High Prevalence Uncertainty], 
    [High Predicted Prevalence \ + Low Prevalence Uncertainty], 
    [High Predicted Prevalence \ + High Prevalence Uncertainty]
  ),

  [High Dose \ + Low Dose Uncertainty], 
    [
      *Interpretation*: Low prevalence in low areas where host biodiversity is expected with high confidence suggests that pathogen control is effective. 
    
      *Action*: Continued collection of prevalence data for population surveillance to ensure that control remains successful over time.
    ],
    [
      *Interpretation*: Confident prediction of high dose with low prevalence indicates control measures are effective, but high uncertainty means more prevalence sampling is necessary. 

      *Action*: Collection of prevalence data to improve the estimate of disease presence in the host community.
    ],
      [
      *Interpretation*: Confident prediction of both high dose and pathogen prevalence indicates high zoonotic hazard, and high risk for spillover if near human populations.

      *Action*: Measures for spillover prevention should be taken, including surveillance in nearby human populations to detect spillover. 
    ],
    [
      *Interpretation*: High prevalence in high dose regions suggests high zoonotic hazard, but high uncertainty in the prevalence indicates the need for more surveillance. 

      *Action*: Sampling for prevalence is necessary to determine if this region truly has the potential for high risk of zoonotic spillover. 
    ],

 [High Dose \ + High Dose Uncertainty], 
    [
      *Interpretation*: Low prevalence indicates control is effective, but the uncertainty in dose emphasizes the need for more data on host composition.  

      *Action*: Targeted host occurrence sampling is necessary to better characterize host distributions and confirm whether the observed low prevalence reflects true absence or data limitations.
    ],
    [
      *Interpretation*: High uncertainty of both dose and prevalence indicates a region of high priority for sampling.

      *Action*: Further collection of both occurrence and prevalence data would reduce uncertainty on both axes.
    ],
    [
      *Interpretation*: Confident predictions of high prevalence among regions with high dose uncertainty indicates the primary uncertainty is whether the pathogen is circulating in many or few hosts.

      *Action*: Both occurrence and prevalence sampling among the diversity of hosts, would improve host community estimates and the understanding of what species the pathogen is circulating in.
    ],
    [
      *Interpretation*: In this case, prevalence is high, but uncertainty in host distributions emphasizes the need for better host composition data.

      *Action*: This scenario calls for prioritizing occurrence data to refine predictions of host presence.
    ],

  [Low Dose \ + Low Dose Uncertainty], 
    [
      *Interpretation*: Given the confident predictions of low dose and prevalence,  further sampling is only valuable for ensuring no change in the composition of hosts over time. 

      *Action*: Continued occurrence sampling to detect biodiversity change is most informative, but should be lower priority than other regions. 
    ],
    [
      *Interpretation*: Confident prediction of low dose, but uncertainty in prevalence indicates there could be active transmission within a small number of hosts.

      *Action*: Further prevalence sampling to determine the presence of pathogens in the sparse host community. 

    ], 
    [
      *Interpretation*: Confident prediction of high prevalence in a known low dose region indicates circulation in one/few host populations.

      *Action*: Further prevalence sampling to determine if the one/few hosts are a reservoir for the pathogen, or undergoing a transient outbreak.
    ],
    [
      *Interpretation*: Uncertainty in prevalence in a region of confident low dose indicates the need to detect if the pathogen is present among one/few hosts. 

      *Action*: Further prevalence sampling to determine if the pathogen is present among the one/few hosts.

    ],

  [Low Dose \ + \ High Dose \ Uncertainty], 

    [
      *Interpretation*: Low prevalence where host biodiversity is uncertain indicates gaps remain in host occurrence data.

      *Action*: Additional species occurrence sampling would confirm host distributions and reduce uncertainty in these regions. 
    ],
    [
      *Interpretation*: High uncertainty in prevalence and host composition indicates the need for both types of data.

      *Action*: Collecting occurrence and to refine host composition first.
    ], 
    [
      *Interpretation*: High prevalence certainty, but uncertainty in dose indicates it is unclear whether the pathogen is present in one/many hosts

      *Action*: Collecting occurrence data to better refine host composition.
    ],
    [
      *Interpretation*: High uncertainty about both dose and prevalence indicates the need for both forms of data
 
      *Action*: This highlights the need to focus sampling efforts on species occurrence data to clarify host composition. 
    ],
)]) <prev-actionable>

