 #single species contributions picrust
#make_unnormalized = get formatted KO abundances for one species
#prmts = get single species PRMT
#contributions - do contributions analysis for all species
##Revision - aggregate to genus level 11/24/2015
##Revision 4/20/2016: single species MUSiCC
#NOTE: need to use OTUs normalized by 16S copy #
#Reorganization: making each of these a function in core_functions
options(stringsAsFactors = F)
library(data.table)
library(reshape2)
library(ggplot2)
library(getopt)
source("core_functions.R")

args = commandArgs(trailingOnly=T)

contrib_file = args[1] #"CD_mice_data/mice_metagenome_contributions.txt" #
data_dir = args[2] #"CD_mice_data/"
results_file = args[3] #"runs_nz4_noMP/mice_AB_out.rda"
out_dir = args[4]
otu_id = args[5] #Usually all
otu_file = args[6] ##Needs to be normalized for 16S copy number!
valueVar = args[7] #RelAbundSample or singleMusicc, can also be determined/reset by previous variables
#data_dir = "BV_data/"
#results_file = "runs_nz4_noMP/bv_val_picrust_out.rda"

cat(paste(otu_id, "\n"))
cat(paste(contrib_file, "\n"))
cat(paste(data_dir, "\n"))

contribs = fread(contrib_file) #To avoid reading in multiple times
all_otus = sort(unique(contribs[,OTU]))
if("make_unnormalized" %in% args){ #Using relative abundance out of all genes
  valueVar = "RelAbundSample"
  contribs = make_unnormalized_single_spec(contribs, otu_id, out_dir)
  if("sum_to_genus" %in% args){
    cat("Summing to genus level, need GG taxonomy file in same directory\n")
    sum_to_genus(contribs, valueVar = valueVar)
  }
  save_contribs_by_species(contribs, valueVar = valueVar, out_dir)
}

if("singleMusicc" %in% args){
  valueVar = "singleMusicc"
  contribs = single_spec_musicc(contribs, otu_file)
  if("sum_to_genus" %in% args){
    cat("Summing to genus level, need GG taxonomy file in same directory\n")
    sum_to_genus(contribs, valueVar = valueVar)
  }
  
  save_contribs_by_species(contribs, valueVar = valueVar, out_dir)

#Comparing this value with MUSiCC values is very close
  #   musicc_genes[KO==ko] 
#   test1 = contribs[Gene==ko,sum(singleMusicc), by=Sample] #Compare this with total. then this value is what we use?
#   test1 = merge(test1, melt(musicc_genes[KO==ko], variable.name="Sample"), by="Sample", all=T)
#   test1[is.na(V1), V1:=0]
#   ggplot(test1, aes(x=value, y=V1, label=Sample)) + geom_text()
}

if("prmts" %in% args){ #Basically just runs get_prmt_scores
  if(otu_id != "all"){
    kos_alone = fread(paste0(data_dir, otu_id, "_KOabund.txt"))
    #kos_alone = fread(paste0(data_dir,otu_id, "_KOabund_musicc.txt"))
    kos_alone = kos_alone[,c(names(kos_alone)[names(kos_alone)!="KO"],"KO"),with=F]
    #load("runs_nz4_noMP/bv_val_picrust_out.rda")
    #load(results_file)
    #Best just to regenerate network
    sub_ko_net = generate_genomic_network(kos_alone[,KO], keggSource = "labKegg",degree_filter = 30)
    prmts_alone = get_prmt_scores(sub_ko_net[[1]], kos_alone)
    write.table(prmts_alone, file = paste0(data_dir, otu_id, "_prmts.txt"), quote=F, row.names=F, sep = "\t")
  }else{
    get_all_singleSpec_prmts(all_otus, valueVar)
  }
}

if(("contributions" %in% args) & otu_id=="all"){
  load(results_file)
  subjects = names(norm_kos)[names(norm_kos)!="KO"]
  prmts_sub_good = get_prmt_scores(ko_net[[1]], norm_kos) # Full community PRMT scores
  prmt_sub_good = prmts_sub_good[node_data[,compound]]
  all_rxns = lapply(node_data[,compound], function(x){ return(ko_net[[3]][Reac==x|Prod==x])})
  all_rxns = lapply(all_rxns,get_non_rev_rxns)
#  if(!all(is.na(match(c("prmts_alone", "all_otus"),ls())))){ #If all of these do not already exist
    load(paste0(out_dir, "all_prmtsAlone_byOTU_",valueVar,".rda"))
    #contribs = fread(contrib_file)
    all_otus = sort(unique(contribs[,OTU]))
 # } 
#for mice abx only - subject id change
#   for(j in 1:length(prmts_alone)){
#     setnames(prmts_alone[[j]], gsub("NoAbx6wkrec","NoAbx6wk", names(prmts_alone[[j]])))
#   }
  spec_contrib=lapply(1:length(node_data[,compound]),prmt_species_contributions_picrust, prmts_sub_good = prmt_sub_good, all_rxns = all_rxns, subjects = subjects, norm_kos = norm_kos, ko_net = ko_net, all_taxa = all_otus, single_spec_prmts = prmts_alone, cor_with=T)
  
  save(spec_contrib, file = paste0(out_dir, "/picrust_spec_contrib_",valueVar,".rda"))  
  
}
