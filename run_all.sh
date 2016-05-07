#!/bin/bash

##Run all taxa-metabolite analyses

HOMEDIR=""

DATADIR1="$HOMEDIR/BV_data/"
DATADIR2="$HOMEDIR/CD_mice_data/"
DATADIR3="$HOMEDIR/Swedish_data/"
DATADIR4="$HOMEDIR/Ecoli_data/"

RUNDIR="$HOMEDIR/runs_out" #output dir for CMP scores

##1 run PICRUSt on Datasets 2 and 3

normalize_by_copy_number.py -f -i "$DATADIR1/Dataset2_otu_table.txt" -o "$DATADIR1/Dataset2_normalized_otus.biom"
predict_metagenomes.py -f -i "$DATADIR1/Dataset2_normalized_otus.biom" -o "$DATADIR1/Dataset2_metagenome_predictions.txt"
metagenome_contributions.py -i "$DATADIR1/Dataset2_normalized_otus.biom" -o "$DATADIR1/bv_val_picrust_metagenome_contributions.txt"

normalize_by_copy_number.py -f -i "$DATADIR2/Dataset3_otu_table.txt" -o "$DATADIR2/Dataset3_normalized_otus.biom"
predict_metagenomes.py -f -i "$DATADIR2/Dataset3_normalized_otus.biom" -o "$DATADIR2/Dataset3_metagenome_predictions.txt"
metagenome_contributions.py -i "$DATADIR2/Dataset3_normalized_otus.biom" -o "$DATADIR2/mice_metagenome_contributions.txt"

#1a clean up picrust output

##2 run in-house version of PICRUSt on Dataset 1
Rscript BV_picrust.R "$DATADIR1/BV_qpcr_good.txt" "$DATADIR1/BV_ref_kos.txt" "$DATADIR1/BV_kos_qpcr_good.txt"


##3 run CMP calculations and comparisons

#BV Dataset 1

gene_file5="$DATADIR1/BV_kos_qpcr_good.txt"
# gene_file5='BV_data/BV_metagenomes_good.txt'
met_file5="$DATADIR1/bv_raw_mets_good.txt"

Rscript runMimosa.R --genefile=$gene_file5 -m $met_file5 -w -p "$RUNDIR/bv_q" -n labKegg -f 30 -z 4 -u 10000


#BV validation Dataset 2
gene_file="$DATADIR1/validation_picrust_good.txt"
met_file="$DATADIR1/validation_mets_good.txt"

run_musicc.py $gene_file -o "$DATADIR1/validation_picrust_musicc.txt" -n -c learn_model -v

Rscript runMimosa.R --genefile="$DATADIR1/validation_picrust_musicc.txt" -m "$DATADIR1/$met_file" -w -p "$RUNDIR/bv_val_picrust" -n labKegg -f 30 -z 4

# #mice ABs Dataset 3

gene_file6="$DATADIR2/mice_picrust_genes_good.txt"
met_file6="$DATADIR2/mice_raw_mets_good.txt"

run_musicc.py "$gene_file6" -o "$DATADIR2/musicc_out.txt" -n -c learn_model -v

Rscript runMimosa.R --genefile="$DATADIR2/musicc_out.txt" -m "$met_file6" -w -p "$RUNDIR/mice_AB" -n labKegg -f 30 -z 4 -u 10000

#Swedish twins Dataset 4
gene_file1="$DATADIR3/genes_good.txt"
met_file1="$DATADIR3/met_data_sub.txt"
met_ids1="$DATADIR3/metabosearch_KEGG_noHMDB.txt"

# #normalize the gene abundance data
run_musicc.py "$gene_file1" -o "$DATADIR3/swedish_musicc_out.txt" -n -c learn_model -v 

Rscript runMimosa.R --genefile="$DATADIR3/swedish_musicc_out.txt" -m "$met_file1"  -w --file_prefix="$RUNDIR/swedish" -n labKegg -f 30 -z 4 -i "$met_ids1"

#E coli dataset
gene_file_ecoli="$HOMEDIR/Ecoli_data/genes_sub_tp14_good.txt"
met_file_ecoli="$HOMEDIR/Ecoli_data/mets_sub_tp14_good.txt"

Rscript runMimosa.R --genefile="$gene_file_ecoli" -m "$met_file_ecoli"  -w --file_prefix="$RUNDIR/ecoli_tp14" -n labKegg -f 30 -z 4


# 4 Get key species contributors for PICRUSt datasets

#Dataset 2
CONTRIBFILE="$DATADIR1/bv_val_picrust_metagenome_contributions.txt"
RESULTSFILE="$RUNDIR/bv_val_picrust_out.rda"

Rscript singleSpec_contributions_picrust.R "$CONTRIBFILE" "$DATADIR1" "$RESULTSFILE" all make_unnormalized prmts contributions save

#Dataset 3
CONTRIBFILE="$DATADIR2/mice_metagenome_contributions.txt"
RESULTSFILE="$RUNDIR/mice_AB_out.rda"

Rscript singleSpec_contributions_picrust.R "$CONTRIBFILE" "$DATADIR2" "$RESULTSFILE" all make_unnormalized prmts contributions save

###Network shuffling tests (submit to cluster)
qsub -t 1-100 $HOMEDIR/shuffle_net.sh

Rscript summarize_net_shuff_tests.R bv_q

##Knit Rmd Document with all core analyses and results
Rscript -e "rmarkdown::render('TaxaMet_paper_figures_all.Rmd')"


