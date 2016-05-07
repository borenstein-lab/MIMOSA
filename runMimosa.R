#7-21-2014
#Functions to create a stoichiometric matrix (EMM) and gene abundance matrix from a metagenomic network, and compare
#with metabolomic data
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
#library(bit64)
options(stringsAsFactors=F)
source("core_functions.R")

spec = matrix(c('genefile','g',1,"character",
              'metfile','m',1,"character",
              'classification','a',0,"logical",
              'write_net','w',0,"logical",
              'file_prefix','p',2,"character",
#              'dir_method','d',2,"character",
              'net_method','n',2,"character",
              'net_file','l',2,"character",
              'degree_filter','f',2,"integer",
              'met_id_file','i',2,"character",
              'minpath_file','k',2,"character",
              'cor_method','c',2,"character",
              'num_permute','u',2,"integer",
              'quant','q',2,"double",
              'nonzero_filt','z',2, "integer"), byrow=T, ncol=4)
  
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
if(!is.null(opt$net_method)) net_method = opt$net_method else net_method = 'load' #load, api, labKegg, or loadNet
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
} else net_file = ""
if(!is.null(opt$num_permute)) num_permute = opt$num_permute else num_permute = 20000
cat(paste("Number of Mantel permutations is ", num_permute,"\n"))
if(!is.null(opt$classification)) {
  runmet2 = T 
  if(!is.null(opt$quant)) quant = opt$quant else quant = 0.5
  cat(paste("Quantile cutoff for classification is"), quant, "\n")
} else runmet2 = F
if(!is.null(opt$nonzero_filt)) nonzero_filt = opt$nonzero_filt else nonzero_filt = 2
cat(paste("Nonzero filter is ", nonzero_filt,"\n"))



if(!runmet2){
  run_all_metabolites(genes, mets, file_prefix = file_prefix, id_met = !is.null(opt$met_id_file), met_id_file = met_id_file, 
                    net_method = net_method, net_file = net_file,
                    correction = "fdr", degree_filter = degree_filter, minpath_file = minpath_file, cor_method = cor_method, nperm = num_permute, nonzero_filter = nonzero_filt)
} else {
  run_all_metabolites2(genes, mets, file_prefix = file_prefix, id_met = !is.null(opt$met_id_file), met_id_file = met_id_file, 
                       net_method = net_method, net_file = net_file,
                       correction = "fdr", degree_filter = degree_filter, minpath_file = minpath_file, quant = quant, nonzero_filter = nonzero_filt)
  
}
