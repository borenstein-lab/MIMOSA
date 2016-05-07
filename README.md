# MIMOSA
Model-based Integration of Metabolite Observations and Species Abundances

A pipeline for joint metabolic model-based analysis of metabolomics measurements and taxonomic composition from microbial communities. 

**Contents**

- core_functions.R: generate community-specific metabolic network model, use to make metabolite predictions and identify consistent and contrasting metabolites

- singleSpec_contributions.R: identify potential taxonomic contributors to metabolite variation

- runMimosa.R: wrapper script for running main analysis

- run_all.sh: example script for re-generating all results from the 4 datasets described in Noecker et al, 2016.

- shuffle_network.R: shuffle the community metabolic network to generate a null model

- network_generation.R: initial generation from community metabolic network template from KEGG database (Kanehisa and Goto, 2000)

**Citation**

Noecker, C., Eng, A., Srinivasan, S., Theriot, C.M., Young, V.B., Jansson, J.K., Fredricks, D.N., and Borenstein, E. (2016). Metabolic Model-Based Integration of Microbiome Taxonomic and Metabolomic Profiles Elucidates Mechanistic Links between Ecological and Metabolic Variation. mSystems 1, e00013–e00015.
