 #single species contributions picrust
#make_unnormalized = get formatted KO abundances for one species
#prmts = get single species PRMT
#contributions - do contributions analysis for all species
##Revision - aggregate to genus level 11/24/2015
##Revision 4/20/2016: single species MUSiCC
#NOTE: need to use OTUs normalized by 16S copy #
#Reorganization: making each of these a function in core_functions
#options(stringsAsFactors = F)
#library(data.table)
#library(reshape2)
#library(ggplot2)
#library(getopt)
#source("core_functions.R")

# args = commandArgs(trailingOnly=T)
#
# contrib_file = args[1] #"CD_mice_data/mice_metagenome_contributions.txt" #
# data_dir = args[2] #"CD_mice_data/"
# results_file = args[3] #"runs_nz4_noMP/mice_AB_out.rda"
# out_dir = args[4]
# otu_id = args[5] #Usually all
# otu_file = args[6] ##Needs to be normalized for 16S copy number!
# valueVar = args[7] #RelAbundSample or singleMusicc, can also be determined/reset by previous variables
# #data_dir = "BV_data/"
# #results_file = "runs_nz4_noMP/bv_val_picrust_out.rda"
#
# cat(paste(otu_id, "\n"))
# cat(paste(contrib_file, "\n"))
# cat(paste(data_dir, "\n"))

#' Evaluate species contributors for a single metabolite with qPCR/species abundance data
#'
#' @import data.table
#' @param j metabolite # (usually from lapply/sapply)
#' @param prmts_sub_good CMP scores for metabolites with abundance data
#' @param all_rxns list of relevant reactions for each metabolite
#' @param subjects vector of subjects
#' @param norm_kos data.table of gene abundances
#' @param ko_net network created by generate_genomic_network
#' @param qpcr original species abundances
#' @param ref_kos gene abundances for each species
#' @param cor_with whether to look at the correlation of CMP scores of each species by itself with the metabolite, or of the whole community with that species removed
#' @return list of potential key species contributors and associated info
#' @examples
#' sapply(1:length(metabolites), prmt_species_contributions, cmp_scores, all_rxns,
#' all_subjects, ko_abunds, ko_net, spec_abunds, ref_kos)
#' @export
prmt_species_contributions = function(j, prmts_sub_good, all_rxns, subjects, norm_kos, ko_net, qpcr, ref_kos, cor_with=F){
  if(!is.null(all_rxns[[j]])){
    if(!cor_with){
      compound = prmts_sub_good[j,compound]
      kos_involved = unique(all_rxns[[j]][Reversible==0,KO])
      #plan - look at correlation between old prmt score and prmt score w/out each species
      species_involved = ref_kos[KO %in% kos_involved, unique(species)]

      species_prmt_cors = sapply(species_involved, function(x){
        spec_without = qpcr[,names(qpcr)!=x, with=F]
        vals_without = kos_from_species(spec_without, ref_kos)
        vals_without = vals_without[KO %in% names(ko_net[[1]])]
        net_mod = ko_net[[1]][,which(names(ko_net[[1]]) %in% vals_without[,KO])]
        prmt_without = data.matrix(net_mod[compound,vals_without[,KO]])%*%data.matrix(vals_without[,subjects,with=F])
        return(cor(as.vector(unlist(prmts_sub_good[compound,subjects,with=F])),as.vector(prmt_without), method="spearman"))
      })
      spec_cors = data.table(Species=species_involved, Cor=species_prmt_cors)
      #save KOs that have a major effect
      species_good = spec_cors[is.na(Cor)|(Cor < 0.5),Species]
      return(list(spec_cors, species_good))
      #  return(ko_cors)
    }else{
      compound = prmts_sub_good[j,compound]
      kos_involved = unique(all_rxns[[j]][Reversible==0,KO])
      species_involved = ref_kos[KO %in% kos_involved, unique(species)]
      species_prmt_cors = sapply(species_involved, function(x){
        spec_alone = qpcr[,names(qpcr) %in% c(x, "Subject","Visit", "Sample", "ID"), with=F]
        vals_alone = kos_from_species(spec_alone, ref_kos)
        vals_alone = vals_alone[KO %in% names(ko_net[[1]])]
        net_mod = ko_net[[1]][,which(names(ko_net[[1]]) %in% vals_alone[,KO])]
        prmt_alone = data.matrix(net_mod[compound,vals_alone[,KO]])%*%data.matrix(vals_alone[,subjects,with=F])
        return(cor(as.vector(unlist(prmts_sub_good[compound,subjects,with=F])),as.vector(prmt_alone), method="spearman"))
      })
      spec_cors = data.table(Species=species_involved, Cor=species_prmt_cors)
      species_good = spec_cors[Cor > 0.5,Species]
      return(list(spec_cors,species_good))
    }
  } else return(NULL)
}

#' Evaluate species contributors for a single metabolite with OTU and PICRUSt data
#'
#' @import data.table
#' @param j metabolite # (usually from lapply/sapply)
#' @param prmts_sub_good CMP scores for metabolites with abundance data
#' @param all_rxns list of relevant reactions for each metabolite
#' @param subjects vector of subjects
#' @param norm_kos data.table of gene abundances
#' @param ko_net network created by generate_genomic_network
#' @param all_taxa vector of OTUs
#' @param single_spec_prmts single-species CMP scores calculated from get_spec_contribs function
#' @param cor_with whether to look at the correlation of CMP scores of each species by itself with the metabolite, or of the whole community with that species removed
#' @return list of 2-item lists for every metabolite - 1st item is data.table of OTUs and correlations, second item is vector of OTUs with correlations > 0.5
#' @examples
#' sapply(1:length(metabolites), prmt_species_contributions_picrust, cmp_scores, all_rxns,
#' all_subjects, ko_abunds, ko_net, spec_abunds, ref_kos)
#' @export
prmt_species_contributions_picrust = function(j, prmts_sub_good, all_rxns, subjects, norm_kos, ko_net, all_taxa, single_spec_prmts, cor_with=F){
  if(!is.null(all_rxns[[j]])){
    compound = prmts_sub_good[j,compound]
    kos_involved = unique(all_rxns[[j]][Reversible==0,KO])
    species_prmt_cors = sapply(1:length(all_taxa), function(x){
      return(cor(as.vector(unlist(prmts_sub_good[compound,subjects,with=F])),as.vector(unlist(single_spec_prmts[[x]][compound,subjects,with=F])), method="spearman"))
    })
    spec_cors = data.table(Species=all_taxa, Cor=species_prmt_cors)
    species_good = spec_cors[Cor > 0.5,Species]
    return(list(spec_cors,species_good))
  } else return(NULL)
}



#' Get contributors for every metabolite and save
#'
#' @param contribs data.table of OTU, genes, samples and abundances
#' @param otu_file File path to OTU abundance matrix
#' @return contribs with abundances normalized to copy number
#' @export
single_spec_musicc = function(contribs, otu_file){
  otus = data.table::fread(otu_file)
  if("#OTU ID" %in% names(otus)) setnames(otus, "#OTU ID", "OTU")
  if("OTU_ID" %in% names(otus)) setnames(otus, "OTU_ID", "OTU")
  otus_rel_abund = otus[,lapply(.SD, function(x){ x/sum(x)})]
  otus_rel_abund[,OTU:=otus[,OTU]]
  otus_rel_abund_melt = data.table::melt(otus_rel_abund, id.var = "OTU", variable.name = "Sample", value.name = "OTURelAbund")
  all_samps = unique(contribs[,Sample])
  #contribs[,RelAbundSample:=CountContributedByOTU/sum(CountContributedByOTU),by=Sample] #relative abundance of this gene from this taxon in all the genes
  contribs = merge(contribs, otus_rel_abund_melt, by=c("OTU", "Sample"), all.x=T)
  all_otus = sort(unique(contribs[,OTU]))
  cat(all_otus)
  all_koAbunds_byOTU = vector("list", length(all_otus))
  all_kos = contribs[,unique(Gene)]
  contribs[,singleMusicc:=OTURelAbund*GeneCountPerGenome] #Relative abundance of OTU times copy number
  return(contribs)
}

#' Get functional relative abundances from a single species' abundances over time
#'
#' @param contribs metagenome_contributions output of PICRUSt
#' @param otu_id specific OTU ID or "all"
#' @param out_dir directory to write to
#' @return gene abundance matrix based on that species only
#' @export
make_unnormalized_single_spec = function(contribs, otu_id = "all", out_dir){
  all_samps = unique(contribs[,Sample])
  contribs[,RelAbundSample:=CountContributedByOTU/sum(CountContributedByOTU),by=Sample] #relative abundance of this gene from this taxon in all the genes
  if(otu_id!="all"){
    contribs = contribs[OTU==otu_id]
    contribs = data.table::dcast.data.table(contribs, Gene~Sample, value.var = "RelAbundSample")
    #contribs = dcast.data.table(contribs, Gene~Sample, value.var = "ContributionPercentOfSample")
    for(j in 1:length(all_samps)){
      if(!(all_samps[j] %in% names(contribs))){
        contribs[,eval(all_samps[j]):=rep(0)]
      }
    }
    setnames(contribs, "Gene","KO")
    contribs = contribs[,c("KO",sort(all_samps)),with=F]
    write.table(contribs, file = paste0(out_dir,otu_id,"_KOabund.txt"), quote=F, row.names=F, sep = "\t")
  }
  return(contribs)
}

#' Sum OTU data to the genus level
#'
#' Needs Greengenes 97_otu_taxonomy.txt to be in the same directory
#'
#' @param contribs metagenome_contributions output of PICRUSt
#' @param valueVar "relAbundSample" or "singleMusicc", abundance metric to use for single taxon gene abundances
#' @return contribs added to the genus level
#' @export
sum_to_genus = function(contribs, valueVar){ #Value var is relAbundSample or singleMusicc
  taxonomy = data.table::fread("97_otu_taxonomy.txt")
  taxonomy[,Genus:=gsub("; s__.*", "", V2)]
  data.table::setnames(taxonomy, c("OTU", "Taxonomy", "Genus"))
  contribs = merge(contribs, taxonomy, by="OTU", all.x=T, all.y=F)
  contribs = contribs[,sum(get(valueVar)), by=list(Gene, Sample, Genus)]
  data.table::setnames(contribs, "Genus", "OTU") #now just treat genera like otus
  data.table::setnames(contribs, "V1", valueVar)
  return(contribs)
}

#' Get contributors for every metabolite and save
#'
#' @param contribs data.table of OTU, genes, samples and abundances
#' @param valueVar "relAbundSample" or "singleMusicc", abundance metric to use for single taxon gene abundances
#' @param out_dir directory to write to
#' @return Null, saves to Rdata file
#' @export
save_contribs_by_species = function(contribs, valueVar, out_dir){
  #Separate contribs by species and save
  all_otus = sort(unique(contribs[,OTU]))
  all_samps = unique(contribs[,as.character(Sample)])
  cat(all_otus)
  all_koAbunds_byOTU = vector("list", length(all_otus))
  for(k in 1:length(all_otus)){
    contribs_sub = contribs[OTU==all_otus[k]]
    contribs_sub = data.table::dcast.data.table(contribs_sub, Gene~Sample, value.var = valueVar)
    #contribs_sub = dcast.data.table(contribs_sub, Gene~Sample, value.var = "ContributionPercentOfSample")
    for(j in 1:length(all_samps)){
      if(!(all_samps[j] %in% names(contribs_sub))){
        contribs_sub[,(as.character(all_samps[j])):=rep(0)]
      }
    }
    data.table::setnames(contribs_sub, "Gene","KO")
    all_koAbunds_byOTU[[k]] = contribs_sub[,c("KO",sort(all_samps)),with=F]
  }
  save(all_koAbunds_byOTU, file = paste0(out_dir,"all_koAbunds_byOTU_",valueVar,".rda"))
}

#' Get single-species CMP scores for every OTU and metabolite
#'
#' @param all_otus vector of OTUs or taxa
#' @param valueVar "relAbundSample" or "singleMusicc", abundance metric to use for single taxon gene abundances
#' @param out_dir Path of directory for loading and saving output
#' @return Null, saves to Rdata file
#' @export
get_all_singleSpec_prmts = function(all_otus, valueVar, out_dir){
  cat(valueVar)
  if(!"all_koAbunds_byOTU" %in% ls()) load(paste0(out_dir, "all_koAbunds_byOTU_",valueVar,".rda"))
  if(length(all_otus) != length(all_koAbunds_byOTU)) stop("Problem!! OTU lists don't match")
  prmts_alone = vector("list", length(all_otus))
  for(k in 1:length(all_otus)){
    sub_ko_net = generate_genomic_network(all_koAbunds_byOTU[[k]][,KO], keggSource = "labKegg", degree_filter = 30)
    prmts_alone[[k]] = get_prmt_scores(sub_ko_net[[1]], all_koAbunds_byOTU[[k]])
  }
  save(prmts_alone, file = paste0(out_dir, "all_prmtsAlone_byOTU_",valueVar,".rda"))
}

#' Perform a series of steps related to identifying single-species contributors.
#'
#' @param contrib_file File where metagenome_contributions PICRUSt output is located
#' @param data_dir Directory where OTU data is located
#' @param results_file File where run_all_metabolites output is located
#' @param out_dir Directory where output will be written too
#' @param otu_id "All" to evaluate all detected OTUs, otherwise, a vector of OTU IDs to test for contribution
#' @param otu_file OTU abundances file
#' @param valueVar type of functional abundances to use for single species calculations - either "singleMusicc" or "RelAbundSample"
#' @param make_unnormalized The next few parameters are binary variables indicating whether to perform each step of the calculation. This is whether to revert to relative abundance values for the gene abundances from single species
#' @param sum_to_genus Whether to sum over OTUs to the genus level (requires Greengenes taxonomy info)
#' @param prmts Whether to calculate single-species CMP scores
#' @param contributions Whether to evaluate contributions based on single-species CMP scores
#' @return No return, saves output to specified file
#' @examples
#' read_files(gene_file, met_file)
#' @export
get_spec_contribs = function(contrib_file, data_dir, results_file, out_dir, otu_id, otu_file, valueVar, make_unnormalized, sum_to_genus, prmts, contributions){
  devtools::load_all()
  contribs = data.table::fread(contrib_file, stringsAsFactors = F) #To avoid reading in multiple times
  all_otus = sort(unique(contribs[,OTU]))
  if(make_unnormalized){ #Using relative abundance out of all genes
    valueVar = "RelAbundSample"
    contribs = make_unnormalized_single_spec(contribs, otu_id, out_dir)
    if("sum_to_genus" %in% args){
      cat("Summing to genus level, need GG taxonomy file in same directory\n")
      sum_to_genus(contribs, valueVar = valueVar)
    }
    save_contribs_by_species(contribs, valueVar = valueVar, out_dir)
  }

  if(valueVar == "singleMusicc"){
    contribs = single_spec_musicc(contribs, otu_file)
    if(sum_to_genus){
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

  if(prmts){ #Basically just runs get_prmt_scores
    if(otu_id != "all"){
      kos_alone = data.table::fread(paste0(data_dir, otu_id, "_KOabund.txt"))
      #kos_alone = data.table::fread(paste0(data_dir,otu_id, "_KOabund_musicc.txt"))
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

  if((contributions == T) & otu_id=="all"){
    load(results_file)
    subjects = names(norm_kos)[names(norm_kos)!="KO"]
    prmts_sub_good = get_prmt_scores(ko_net[[1]], norm_kos) # Full community PRMT scores
    prmt_sub_good = prmts_sub_good[node_data[,compound]]
    all_rxns = lapply(node_data[,compound], function(x){ return(ko_net[[3]][Reac==x|Prod==x])})
    all_rxns = lapply(all_rxns,get_non_rev_rxns)
    #  if(!all(is.na(match(c("prmts_alone", "all_otus"),ls())))){ #If all of these do not already exist
    load(paste0(out_dir, "all_prmtsAlone_byOTU_",valueVar,".rda"))
    #contribs = data.table::fread(contrib_file)
    all_otus = sort(unique(contribs[,OTU]))
    # }
    #for mice abx only - subject id change
    #   for(j in 1:length(prmts_alone)){
    #     setnames(prmts_alone[[j]], gsub("NoAbx6wkrec","NoAbx6wk", names(prmts_alone[[j]])))
    #   }
    spec_contrib=lapply(1:length(node_data[,compound]),prmt_species_contributions_picrust, prmts_sub_good = prmt_sub_good, all_rxns = all_rxns, subjects = subjects, norm_kos = norm_kos, ko_net = ko_net, all_taxa = all_otus, single_spec_prmts = prmts_alone, cor_with=T)

    save(spec_contrib, file = paste0(out_dir, "/picrust_spec_contrib_",valueVar,".rda"))
  }
  return(NULL)
}

