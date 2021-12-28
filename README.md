Arctic host-parasite dynamics under climate change more influenced by
behaviour than thermal performance
================
Stephanie Peacock
December 28, 2021

##### “OsterBou”= Ostertagia + cariBou

##### “pop” = population dynamics

This repository contains data and code to reproduce the results
presented in the submission by Peacock et al. to Royal Society Open
Science, titled “Arctic host-parasite dynamics under climate change more
influenced by behaviour than thermal performance”. The following folders
contain the primary code and data presented:

1.  `MTE_dataAnalysis` contains data (`Ostertagia_data.csv`) from lab
    experiments to estimate the temperature sensitivity of development
    and mortality of free-living larvae from the parasitic nematode,
    Ostertagia gruehneri, and code to estimate the parameters in thermal
    response curves for both development and mortality from those data
    (`MTE_fitting.R`)

2.  `simulations` contains code for simulating a migratory host-parasite
    model (Peacock et al. 2018, 2020) tailored to barren-ground caribou
    and Ostertagia, including the thermal response curves described in
    (1). The file `PaperSims.R` contains code to run the simulations and
    save output as .rds files, while `PaperFigs.R` reads those .rds
    files to produce the figures included in the manuscript.

Other folders containing information that is peripheral to the work
submitted for publication include:

-   `parameterization` contains background information, including some
    code and figures, exploring functions for parameters in the
    simulations other than the free-living larvae mortality and
    development described in (1) (Appendix B).

-   `tempData` includes data and code used to construct the ground
    surface temperature profiles along the caribou migration route as
    well as the projected temperature anomalies under climate change
    (Appendix C).

## Contact

Questions about this code can be directed to Steph Peacock
<stephanie.j.peacock@gmail.com>
