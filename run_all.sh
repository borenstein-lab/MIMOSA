#!/bin/bash

## Run all taxa-metabolite analyses
## The example OTU table and metabolite measurements used here can be downloaded from elbo.gs.washington.edu/download.html

HOMEDIR=""

DATADIR1=""

RUNDIR="$HOMEDIR/runs_out" #output dir for CMP scores

##1 run PICRUSt

normalize_by_copy_number.py -f -i "$DATADIR1/Dataset2_otu_table.txt" -o "$DATADIR1/Dataset2_normalized_otus.biom"
predict_metagenomes.py -f -i "$DATADIR1/Dataset2_normalized_otus.biom" -o "$DATADIR1/Dataset2_metagenome_predictions.txt"
metagenome_contributions.py -i "$DATADIR1/Dataset2_normalized_otus.biom" -o "$DATADIR1/Dataset2_metagenome_contributions.txt"

##2 run CMP calculations and comparisons

gene_file="$DATADIR1/Dataset2_metagenome_predictions.txt"
met_file="$DATADIR1/Dataset2_mets.txt"

run_musicc.py "$gene_file" -o "$DATADIR1/Dataset2_picrust_musicc.txt" -n -c learn_model -v

#Using KEGGREST to build KEGG community network, can be slow
Rscript runMimosa.R --genefile="$DATADIR1/Dataset2_picrust_musicc.txt" -m "$met_file" -w -p "$RUNDIR/Dataset2_bv" -n KeggTemplate -f 30 -z 4 -e "$HOMEDIR/reaction_mapformula.lst"

#OR Using downloaded KEGG files to build community network
Rscript runMimosa.R --genefile="$DATADIR1/Dataset2_picrust_musicc.txt" -m "$met_file" -w -p "$RUNDIR/Dataset2_bv" -n KeggTemplate -f 30 -z 4 -e "$HOMEDIR/reaction_mapformula.lst" -r "ko_reaction.list" -x "reaction"


### Options for runMimosa.R:
# -g,--genefile gene abundance file path
# -m,--metfile metabolite abundance file path
# -a,--classification whether to compare with a threshold instead of making all pairwise comparisons with a Mantel test
# -c,--cor_method Whether to use Spearman or Pearson correlations for the Mantel test
# -u,--num_permute How many permutations to perform for the Mantel test
# -q,--quant If classification is selected, what quantile to use as the threshold for "high" vs "low" metabolite levels
# -w,--write_net whether to save the dataset-specific community metabolic network model as part of the output
# -p,--file_prefix file path prefix for saving output files
# -n,--net_method method to generate community metabolic network model: either loadNet or KeggTemplate
# -f,--degree_filter Connectivity threshold, metabolites connected by the model to this many genes or more will be filtered from the analysis
# -z,--nonzero_filt Nonzero filter, metabolites with fewer than this number of samples with nonzero metabolite concentrations or metabolic potential scores will be filtered from the analysis
# -i,--met_id_file File of putative KEGG compound IDs associated with metabolite peak masses for unidentified metabolomics data, output of MetaboSearch
# -e,--mapformula_file File path to reaction_mapformula.lst file, containing information on reactions annotated in KEGG pathways
# -r,--ko_rxn_file File path to links between KO gene families and reactions. Optional - if not provided will use the KEGGREST package to obtain this information from the KEGG API
# -x,rxn_annots_file File path to KEGG Reaction annotations. Optional - if not provided will use the KEGGREST package to obtain this information from the KEGG API
# -b,--keggFile Rda file path of pregenerated network template, containing "all_kegg" object
# -o,--contribs_file file path of PICRUSt-style KO abundances from each OTU, output of metagenome_contributions.py
