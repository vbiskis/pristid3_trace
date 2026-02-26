# pristid3_trace
Elemental analysis of sawfish teeth and vertebrae

------------------------------------------------------------------------

This repo contains R code re: the analyses for the manuscript: 

> Biskis, V. N., Reis-Santos, P., Holmes, B. J., Gillanders, B. M., Townsend, K. A. & Wueringer, B. E. (*2026*). [The potential of elemental analysis in determining origin of sawfish rostra.] In Prep.

Additional output is outlined in the completed PhD thesis:

> Biskis, V. N. (2026). [Ch. 4: The potential of elemental analysis...] In: Threatened sawfishes of Australia: Spatial ecology and species composition through time. PhD Dissertation. *University of the Sunshine Coast*, (https://doi.org/10.25907/01001)

------------------------------------------------------------------------

## Author: Nikki Biskis

PhD UniSC | SARA
Date: 2025-12-10

------------------------------------------------------------------------

## Contents

pristid2_maps

    └── scripts             <- code used in the analysis 
    └── helpers             <- source code written by author VNB   
        ├── functions       <- stats functions written by author VNB, unless otherwise stated 
    └── figs                <- .png files of figures in main text   
        ├── supps           <- additional visuals for the manuscript appendix
        ├── thesis_vs       <- figures that appear in the thesis, but not the published journal article

------------------------------------------------------------------------

## Notes

Code was initially written in August - October 2024, with analysis continuing into November 2025. Edits prior to thesis submission in June 2025.

In February, 2026, the following modifications were done clarity/reproducibility:

1. Functions created for repetitive (multi species or regions) analyses
     - This is noted at the beginning of each script
    
2. Where a large object created (ggplot), saved as rds for later loading

3. Further developed DTW methodology for publication.
 
## Questions & feedback?

Please submit an issue, or email your questions to Nikki Biskis: nikkibiskis@gmail.com

---
