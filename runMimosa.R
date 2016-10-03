#Create a stoichiometric matrix (EMM) and gene abundance matrix from a metagenomic network, and compare with metabolomic data
#assumes already-processed KO and metabolite abundance datasets

library(KEGGREST)
library(data.table)
library(RColorBrewer)
library(vegan)
library(stringr)
library(qvalue)
library(getopt)
library(ggplot2)
library(reshape2)
library(mimosa)
options(stringsAsFactors=F, curl_interrupt = F)
#source("core_functions.R")

spec = matrix(c('genefile','g',1,"character",
              'metfile','m',1,"character",
              'classification','a',0,"logical",
              'write_net','w',0,"logical",
              'file_prefix','p',2,"character",
              'net_method','n',2,"character",
              'net_file','l',2,"character",
              'degree_filter','f',2,"integer",
              'met_id_file','i',2,"character",
              'minpath_file','k',2,"character",
              'cor_method','c',2,"character",
              'num_permute','u',2,"integer",
              'quant','q',2,"double",
              'nonzero_filt','z',2, "integer",
              'mapformula_file','e',2, "character",
              'ko_rxn_file','r',2, "character",
              'rxn_annots_file','x', 2, "character",
              'contribs_file', 'o', 2, "character",
              'keggFile', 'b', 2, "character"), byrow=T, ncol=4)

opt = getopt(spec, opt = commandArgs(TRUE))
datasets = read_files(opt$genefile, opt$metfile)
genes = datasets[[1]]
mets = datasets[[2]]
if(opt$write_net) write_net = T else write_net = F
cat(paste("Write_net is ",write_net,"\n"))
if(!is.null(opt$file_prefix)) file_prefix = opt$file_prefix else file_prefix = 'net1'
cat(paste("File prefix is ", file_prefix,"\n"))
#if(!is.null(opt$dir_method)) dir_method = opt$dir_method else dir_method = 'standard' #standard or minpath
#cat(paste("Dir method is ", dir_method,"\n"))
if(!is.null(opt$net_method)) net_method = opt$net_method else net_method = 'load' #loadNet or KeggTemplate
cat(paste("Net method is ", net_method,"\n"))
if(!is.null(opt$degree_filter)) degree_filter = opt$degree_filter else degree_filter = 0
cat(paste("Degree filter is", degree_filter,"\n"))
if(!is.null(opt$met_id_file)) met_id_file = opt$met_id_file else met_id_file = ''
cat(paste("Met id file is ", met_id_file,"\n"))
if(!is.null(opt$minpath_file)) minpath_file = opt$minpath_file else minpath_file = ''
cat(paste("Minpath file is ", minpath_file,"\n"))
if(!is.null(opt$cor_method)) cor_method = opt$cor_method else cor_method = "spearman"
cat(paste("Mantel correlation method is ",cor_method,"\n"))
if(net_method == "loadNet"){
  if(!is.null(opt$net_file)) net_file = opt$net_file
  cat(paste("Network file is ",net_file,"\n"))
} else if(net_method == "KeggTemplate"){
  if(is.null(opt$mapformula_file)){
    stop("Need mapformula file!")
  } else{
    #Set up generic network info
    cat("Generating network template from KEGG...\n")
    if(is.null(opt$ko_rxn_file) & is.null(opt$keggFile)){ ##Use KEGGREST
      all_kegg = get_kegg_reaction_info("KEGGREST", kolist = genes[,KO])
    } else{
      if(!is.null(opt$keggFile)){
        load(opt$keggFile) #should include correctly formatted object named "all_kegg"
      } else if(is.null(opt$rxn_annots_file)){
        stop("If using downloaded KEGG files for reaction info, must provide a file linking KOs to reactions and a file with reaction annotations")
      } else{
        all_kegg = get_kegg_reaction_info(opt$ko_rxn_file, opt$rxn_annots_file)
      }
    }
    rxn_table = generate_network_template_kegg(opt$mapformula_file, all_kegg)
    cat("Got community network!\n")
  }
}
if(!is.null(opt$num_permute)) num_permute = opt$num_permute else num_permute = 20000
cat(paste("Number of Mantel permutations is ", num_permute,"\n"))
if(!is.null(opt$classification)) {
  runmet2 = T
  if(!is.null(opt$quant)) quant = opt$quant else quant = 0.5
  cat(paste("Quantile cutoff for classification is"), quant, "\n")
} else runmet2 = F
if(!is.null(opt$nonzero_filt)) nonzero_filt = opt$nonzero_filt else nonzero_filt = 3
cat(paste("Nonzero filter is ", nonzero_filt,"\n"))

cat("Running main MIMOSA analysis\n")
if(!runmet2){
  run_all_metabolites(genes, mets, file_prefix = file_prefix, id_met = !is.null(opt$met_id_file), met_id_file = met_id_file,
                    net_method = net_method, net_file = net_file, rxn_table_source = rxn_table,
                    correction = "fdr", degree_filter = degree_filter, minpath_file = minpath_file, cor_method = cor_method, nperm = num_permute, nonzero_filter = nonzero_filt)
} else {
  run_all_metabolites2(genes, mets, file_prefix = file_prefix, id_met = !is.null(opt$met_id_file), met_id_file = met_id_file,
                       net_method = net_method, net_file = net_file,  rxn_table_source = rxn_table,
                       correction = "fdr", degree_filter = degree_filter, minpath_file = minpath_file, quant = quant, nonzero_filter = nonzero_filt)

}

### Get potential key species contributors, assumes PICRUSt was used to generate metagenome predictions
if(!is.null(opt$contribs_file)){
  cat("Getting potential species contributors to metabolite variation\n")
  get_spec_contribs(opt$contribs_file, data_dir = getwd(), results_file = paste0(file_prefix, "_out.rda"), out_dir = getwd(), otu_file = "Dataset2_otu_table.txt", otu_id = "all", valueVar = "singleMusicc", make_unnormalized = F, sum_to_genus = F, cmps = T, contributions = T)
}

### Get key gene/reaction contributors across all species
cat("Getting potential gene and reaction contributors to metabolite variation\n")
load(paste0(file_prefix, "_out.rda"))
good_mets = node_data[,compound]
subjects = names(mets)[names(mets) != "KEGG"]
cmps = get_cmp_scores(ko_net[[1]], norm_kos)
cmps_sub_good = cmps[compound %in% good_mets]
all_rxns = lapply(good_mets, function(x){ return(get_non_rev_rxns(ko_net[[3]][Reac==x | Prod==x]))})
all_gene_contribs = lapply(1:length(good_mets), gene_contributions, cmps_sub_good = cmps_sub_good, all_rxns = all_rxns,
       subjects=subjects, norm_kos = norm_kos, ko_net = ko_net)
save(all_gene_contribs, file = paste0(file_prefix, "_geneContribs.rda"))

### Run 500 network shuffling tests
for(j in 1:500){
  cat(paste0("Running network shuffling test, iteration ", j ,"\n"))
  run_shuffle(paste0(file_prefix, "_out.rda"), id_num = j)
}
