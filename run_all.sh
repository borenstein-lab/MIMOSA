#!/bin/bash

## Run all taxa-metabolite analyses
## The example OTU table and metabolite measurements used here can be downloaded from elbo.gs.washington.edu/download.html

HOMEDIR="" #location of KEGG data

DATADIR1="" #location of metabolite data

RUNDIR="" #output dir for CMP scores

##1 run PICRUSt

normalize_by_copy_number.py -f -i "$DATADIR1/Dataset2_otu_table.txt" -o "$DATADIR1/Dataset2_normalized_otus.biom"
predict_metagenomes.py -f -i "$DATADIR1/Dataset2_normalized_otus.biom" -o "$DATADIR1/Dataset2_metagenome_predictions.txt"
metagenome_contributions.py -i "$DATADIR1/Dataset2_normalized_otus.biom" -o "$DATADIR1/Dataset2_metagenome_contributions.txt"

###Remove extra column of annotations
rev Dataset2_metagenome_predictions.txt | cut -d "  " -f 2- | rev > tmp
mv tmp Dataset2_metagenome_predictions.txt

##2 normalize gene abundances with MUSiCC

gene_file="$DATADIR1/Dataset2_metagenome_predictions.txt"
met_file="$DATADIR1/Dataset2_mets.txt"

run_musicc.py "$gene_file" -o "$DATADIR1/Dataset2_picrust_musicc.txt" -n -v

#Using KEGGREST to build KEGG community network, can be slow
Rscript runMimosa.R --genefile="Dataset2_picrust_musicc.txt" -m "Dataset2_mets.txt" --contribs_file="Dataset2_metagenome_contributions.txt" --file_prefix="Dataset2_bv" --mapformula_file="keggPath/reaction_mapformula.lst" --nonzero_filt=4 

#OR Using downloaded KEGG files to build community network
Rscript runMimosa.R --genefile="Dataset2_picrust_musicc.txt" --metfile="Dataset2_mets.txt" --contribs_file="Dataset2_metagenome_contributions.txt" --file_prefix="Dataset2_bv" --mapformula_file="keggPath/reaction_mapformula.lst" --ko_rxn_file="keggPath/ko_reaction.list" --rxn_annots_file="keggPath/reaction" --nonzero_filt=4 


### Options for runMimosa.R:
# -g,--genefile gene abundance file path
# -m,--metfile metabolite abundance file path
# -o,--contribs_file file path of PICRUSt-style KO abundances from each OTU, output of metagenome_contributions.py
# -e,--mapformula_file File path to reaction_mapformula.lst file, containing information on reactions annotated in KEGG pathways
# -p,--file_prefix file path prefix for saving output files
# -r,--ko_rxn_file File path to links between KO gene families and reactions. Optional - if not provided will use the KEGGREST package to obtain this information from the KEGG API
# -x,rxn_annots_file File path to KEGG Reaction annotations. Optional - if not provided will use the KEGGREST package to obtain this information from the KEGG API
# -b,--keggFile Rda file path of pregenerated network template, containing "all_kegg" object
# -c,--cor_method Whether to use Spearman or Pearson correlations for the Mantel test
# -u,--num_permute How many permutations to perform for the Mantel test
# -f,--degree_filter Connectivity threshold, metabolites connected by the model to this many genes or more will be filtered from the analysis
# -z,--nonzero_filt Nonzero filter, metabolites with fewer than this number of samples with nonzero metabolite concentrations or metabolic potential scores will be filtered from the analysis
# -a,--classification whether to compare with a threshold instead of making all pairwise comparisons with a Mantel test
# -q,--quant If classification is selected, what quantile to use as the threshold for "high" vs "low" metabolite levels
# -i,--met_id_file File of putative KEGG compound IDs associated with metabolite peak masses for unidentified metabolomics data, output of MetaboSearch
# -w,--write_net whether to save the dataset-specific community metabolic network model as part of the output
# -n,--net_method method to generate community metabolic network model: either loadNet or KeggTemplate
