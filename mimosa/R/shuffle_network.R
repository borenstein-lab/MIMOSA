##Network shuffling code
#library(data.table)
#library(reshape2)
#library(vegan)
#library(qvalue)
#options(stringsAsFactors = FALSE)

<<<<<<< HEAD
=======
>>>>>>> origin/master
#source(paste0(homedir, "PRMT_functions_clean.R"))

##Function to randomize a network
##randomly sample 2 edges and if it works, switch products
randomize_net=function(netw, n_reps){
  revnet = get_non_rev_rxns(netw, all_rxns=F)
  network2=revnet #we'll expand the reversible edges back out later
  #Don't want the switching to be biased towards reversible edges
  n_edges = dim(revnet)[1]
  m = 1
  while(m < n_reps){
    rand_is=sample(n_edges,size=2)
    rand_edges=network2[rand_is]
    if(all(rand_edges[,Reversible]==0) | all(rand_edges[,Reversible]==1)){ ##if both are single edges or both are reversible
      ##check that switched connection does not already exist
      comps=unique(c(rand_edges[,Prod], rand_edges[,Reac]))
      check=network2[Reac %in% comps & Prod %in% comps & KO %in% rand_edges[,KO]]
      if(dim(check)[1]==dim(rand_edges)[1] & length(comps)==4){
        #make the switch
        ind1 = which(network2[,KO]==rand_edges[1,KO] & network2[,Reac]==rand_edges[1,Reac] & network2[,Prod]==rand_edges[1,Prod])
        network2[ind1, Prod:=rand_edges[2,Prod]]
        network2[ind1, stoichProd:=rand_edges[2,stoichProd]]
        ind2 = which(network2[,KO]==rand_edges[2,KO] & network2[,Reac]==rand_edges[2,Reac] & network2[,Prod]==rand_edges[2,Prod])
        network2[ind2, Prod:=rand_edges[1,Prod]]
        network2[ind2, stoichProd:=rand_edges[1,stoichProd]]
        m=m+1
      }
    }
  }
  #Add back reversible complements
  for(j in 1:length(network2[,KO])){
    if(network2[j,Reversible]==1){
      network2 = rbind(network2, data.table(KO=network2[j,KO], Reac=network2[j,Prod], Prod=network2[j,Reac], stoichReac=network2[j,stoichProd], stoichProd=network2[j,stoichReac], Reversible=1))
    }
  }
  network2[,stoichReac:=as.numeric(stoichReac)]
  network2[,stoichProd:=as.numeric(stoichProd)]
  return(network2)
}

make_network_matrix = function(rxn_table){
  cmpds = unique(c(rxn_table[,Reac], rxn_table[,Prod]))
  goodkos = unique(rxn_table[,KO])
  network_mat = matrix(rep(0), nrow = length(cmpds), ncol = length(goodkos))
  stoich_mat = matrix(rep(NA), nrow = length(cmpds), ncol = length(goodkos))
  for(j in 1:length(rxn_table[,KO])){
    foo1 = match(rxn_table[j,Reac], cmpds)
    foo2 = match(rxn_table[j,Prod], cmpds)
    fooko = match(rxn_table[j,KO], goodkos)
    network_mat[foo1, fooko] = network_mat[foo1, fooko] - rxn_table[j,stoichReac]
    network_mat[foo2, fooko] = network_mat[foo2, fooko] + rxn_table[j, stoichProd]
    if(is.na(stoich_mat[foo1,fooko])) stoich_mat[foo1,fooko] = -1*rxn_table[j,stoichReac] else stoich_mat[foo1,fooko] = stoich_mat[foo1, fooko] - rxn_table[j,stoichReac]
    if(is.na(stoich_mat[foo2,fooko])) stoich_mat[foo2,fooko] = 1*rxn_table[j,stoichProd] else stoich_mat[foo2,fooko] = stoich_mat[foo2, fooko] + rxn_table[j,stoichProd]
  }
  negsums = apply(network_mat, 1, function(x){ abs(sum(x[x < 0]))})
  possums = apply(network_mat, 1, function(x){ sum(x[x > 0])})
  for(j in 1:length(cmpds)){
    negkos = which(network_mat[j,] < 0)
    poskos = which(network_mat[j,] > 0)
    if(length(negkos) > 0) network_mat[j,negkos] = network_mat[j,][negkos]/negsums[j]
    if(length(poskos) > 0) network_mat[j,poskos] = network_mat[j,][poskos]/possums[j]
  }
  network_mat = data.frame(network_mat)
  names(network_mat) = goodkos
  row.names(network_mat) = cmpds
  stoich_mat = data.frame(stoich_mat)
  names(stoich_mat) = goodkos
  row.names(stoich_mat) = cmpds
  return(list(network_mat, stoich_mat))
}


# args = commandArgs(trailingOnly = T)
# out_file = args[1]
# #out_file = "runs_cleanNet_noMP/bv_q_out.rda"
# load(out_file)
#
# nonzero_filter = 2
#
# id_num = args[2]
#
# rxn_table = ko_net[[3]]
# n_iter = 5000
# random_net = randomize_net(rxn_table, n_iter)
# net_mats = make_network_matrix(random_net)
#
# prmts = get_prmt_scores(net_mats[[1]], norm_kos)
#
# metIDs = mets[,KEGG]
# shared_mets = metIDs[metIDs %in% row.names(net_mats[[1]])]
# mets_shared_only = mets[shared_mets]
# prmts_shared_only = prmts[shared_mets]
#
# good_data = which(apply(prmts_shared_only, 1, function(x){ length(x[as.numeric(x)!=0]) > nonzero_filter + 1}) & apply(mets_shared_only, 1, function(x){ length(x[as.numeric(x)!=0]) > nonzero_filter + 1}))
# mets_shared_only = mets_shared_only[good_data]
# prmts_shared_only = prmts_shared_only[good_data]
# shared_mets = shared_mets[good_data]
#
# compare_pos = sapply(1:length(shared_mets), function(x){
#   compare_met(shared_mets[x], shared_mets[x], mets_shared_only, prmts_shared_only, posneg = "pos", nperm = 10000)
# })
# compare_neg = sapply(1:length(shared_mets), function(x){
#   compare_met(shared_mets[x], shared_mets[x], mets_shared_only, prmts_shared_only, posneg = "neg", nperm = 10000)
# })
#
# corrected_pos = correct(compare_pos, method="fdr")
# corrected_neg = correct(compare_neg, method = "fdr")
# total_pos1 = length(corrected_pos[corrected_pos < 0.1 & compare_pos < 0.05 & !is.na(corrected_pos)])
# total_neg1 = length(corrected_neg[corrected_neg < 0.1 & compare_neg < 0.05 & !is.na(corrected_neg)])
# total_pos2 = length(corrected_pos[corrected_pos < 0.01 & compare_pos < 0.01 & !is.na(corrected_pos)])
# total_neg2 = length(corrected_neg[corrected_neg < 0.01 & compare_neg < 0.01 & !is.na(corrected_neg)])
#
# met_table = data.table(compound=shared_mets, Pos=corrected_pos, Neg=corrected_neg, Iter=rep(id_num, length(shared_mets)))
#
# file_prefix = paste0(gsub("_out.rda", "", out_file), id_num)
# cat(file_prefix)
# write.table(data.frame(total_pos1, total_neg1, total_pos2, total_neg2), file = paste0(file_prefix, "_shuff_network.txt"), sep = "\t", quote=F, row.names=F, col.names=F)
#
# write.table(met_table, file = paste0(file_prefix,"_shuff_network_mets.txt"), sep = "\t", quote=F, row.names=F, col.names=F)
#



