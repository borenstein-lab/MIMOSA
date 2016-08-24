#PRMT_functions.R
#helper script defining functions needed for PRMT analysis

#Add-ons to ggplot2
# guides_merge <- function(gdefs) {
#   gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
#   tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
# }
# environment(guides_merge) <- environment(ggplot)
# assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")
#
# g_legend<-function(a.gplot){
#   tmp <- ggplot_gtable(ggplot_build(a.gplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   return(legend)}
#
# g_title = function(b.gplot){
#   tmp <- ggplot_gtable(ggplot_build(b.gplot))
#   title_i <- which(sapply(tmp$grobs, function(x){grepl("plot.title", x$name)}))
#   title <- tmp$grobs[[title_i]]
#   return(title)
# }
#
# ggsave <- ggplot2::ggsave; body(ggsave) <- body(ggplot2::ggsave)[-2]

#' Read in processed files and assign keys to datasets. Sets NA metabolite values to 0.
#'
#' @import data.table
#' @param genefile File where gene data is located.
#' @param metfile File where metaboltie data is located.
#' @return A list in which the first item is a data.table of gene abundances and the second is a data.table of metabolite abundances.
#' @examples
#' read_files(gene_file, met_file)
#' @export
read_files = function(genefile, metfile){
  genes = fread(genefile, header=T, sep="\t")
  setkey(genes,KO)
  mets = fread(metfile, header=T, sep="\t")
  if("KEGG" %in% names(mets)){
	 mets = mets[,c(subjects,"KEGG"), with=F]
  } else mets = mets[,c("Mass", subjects), with=F]
   #Set NAs to 0
  for(j in names(mets)){
    set(mets,which(is.na(mets[[j]])),j,0)
  }
  if("KEGG" %in% names(mets)) setkey(mets,KEGG) #2 possibilities for metabolite file format
  #save only samples that have both kinds of data and put datasets in the same order of subjects/samples
  subjects = sort(intersect(names(genes), names(mets)))
  genes = genes[,c("KO", subjects), with=F]
  if("KEGG" %in% names(mets)) mets = mets[,c(subjects,"KEGG"), with=F]
  else mets = mets[,c("Mass", subjects), with=F]
  return(list(genes, mets))
}

#' Create a community metabolic network model using a few different methods.
#'
#' @import data.table
#' @param kos Genes to include in network (KEGG Orthology IDs)
#' @param keggSource source of network information, currently can be "labKEGG", "loadNet", or "metacyc"
#' @param degree_filter Compounds connected to this number of KOs or more will be filtered from the network
#' @param minpath_file file of minimal pathways to use for core network
#' @param normalize Whether to normalize output matrix to show relative impacts of each gene on each compound
#' @param networkFile If keggSource is "loadNet", file that template network should be loaded from
#' @return List containing 3 different versions of the same network: an adjacency matrix, an adjacency matrix that differentiates between genes with a neutral effect on a compound and no effect (0 vs NA), and an edge list.
#' @examples
#' generate_genomic_network(kos, "labKegg", degree_filter = 30)
#' @export
generate_genomic_network = function(kos, keggSource = "labKegg", degree_filter = 0,
                                    minpath_file = '', normalize = T, networkFile=""){
  #keggSource must be either "labKegg", or "loadNet" for if the network has already been generated, or "metacyc"
  #degree filter says whether to filter hub compounds with a degree greater than this value
  if(keggSource == "loadNet"){
    load(networkFile)
    return(allnet)
  } else if(keggSource == "labKegg") { #load sharon's network and just grab subset
    #cat("Using mapformula/labKegg\n")
    rxn_table = fread("ko_rxn_map_all_info_filtered.txt", colClasses = c(rep("character",6), rep("numeric",2)))
    rxn_table = rxn_table[KO %in% kos]

    #rxn_table[,rxn_id:=rxn_ids2]
    if(minpath_file!=''){
      minpaths = fread(minpath_file, colClasses="character")
      setnames(minpaths,"Path")
      #for reactions in the minpath set we are only saving the info from those minimal pathways,
      #but we are also saving the reactions not in the minpath set
      rxn_table2 = rxn_table[Path %in% minpaths[,Path]]
      rxn_table_extra = rxn_table[!(Rxn %in% rxn_table2[,Rxn])]
      rxn_table = rbind(rxn_table2, rxn_table_extra)
    }
    #No longer need extra info
    rxn_table[,Path:=NULL]
    rxn_table[,ReacProd:=NULL]
    rxn_table[,Rxn:=NULL]
    rxn_table = unique(rxn_table)
    if(degree_filter != 0){
      #cat("Filtering currency metabolites\n")
      cmpds = unique(c(rxn_table[,Prod], rxn_table[,Reac]))
      degree = sapply(cmpds, function(x){
        rxn_table[Prod==x | Reac==x, length(unique(KO))]
      })
      cmpds = cmpds[degree < degree_filter]
      degree = degree[degree < degree_filter]
      rxn_table = rxn_table[(Prod %in% cmpds) & (Reac %in% cmpds)]
      #redefine cmpds
      cmpds = unique(c(rxn_table[,Prod], rxn_table[,Reac]))
    }
    goodkos = unique(rxn_table[,KO])
    #expand into adjacency and stoichiometry matrices
    #cat("Making stoichiometric matrix\n")
    network_mat = matrix(rep(0), nrow = length(cmpds), ncol = length(goodkos))
    stoich_mat = matrix(rep(NA), nrow = length(cmpds), ncol = length(goodkos))
    for(j in 1:length(rxn_table[,KO])){
      foo1 = match(rxn_table[j,Reac], cmpds)
      foo2 = match(rxn_table[j,Prod], cmpds)
      fooko = match(rxn_table[j,KO], goodkos)
      network_mat[foo1, fooko] = network_mat[foo1, fooko] - rxn_table[j,stoichReac]
      network_mat[foo2, fooko] = network_mat[foo2, fooko] + rxn_table[j, stoichProd]
      if(is.na(stoich_mat[foo1,fooko])) stoich_mat[foo1,fooko] = -1*rxn_table[j,stoichReac] else stoich_mat[foo1,fooko] = stoich_mat[foo1, fooko] - rxn_table[j,stoichReac]
      if(is.na(stoich_mat[foo2,fooko])) stoich_mat[foo2,fooko] = 1*rxn_table[j,stoichProd] else stoich_mat[foo2,fooko] = stoich_mat[foo2, fooko] + rxn_table[j,stoichProd] ##Luckily this doesn't affect anything
    }
      negsums = apply(network_mat, 1, function(x){ sum(x[x < 0])})
      possums = apply(network_mat, 1, function(x){ sum(x[x > 0])})
      metlen = dim(network_mat)[1]
      for(j in 1:metlen){
        negkos = which(network_mat[j,] < 0)
        poskos = which(network_mat[j,] > 0)
        if(length(negkos) > 0) network_mat[j,negkos] = network_mat[j,][negkos]/abs(negsums[j])
        if(length(poskos) > 0) network_mat[j,poskos] = network_mat[j,][poskos]/possums[j]
      }
    network_mat = data.frame(network_mat)
    names(network_mat) = goodkos
    row.names(network_mat) = cmpds
    stoich_mat = data.frame(stoich_mat)
    names(stoich_mat) = goodkos
    row.names(stoich_mat) = cmpds
    #cat("Done with network!\n")
    return(list(network_mat, stoich_mat, rxn_table))
  } else if(keggSource == "metacyc"){
    #Get KEGG reaction IDs for these KOs
    load("KeggReactions.rda") # should give all rxns
    goodrxns = sapply(all_kegg$Reaction_info, function(x){ any(names(x$ORTHOLOGY) %in% kos)})
    rxns = all_kegg$Reactions[goodrxns]
    rxn_ids = rxns
    rxn_info = all_kegg$Reaction_info[goodrxns]

    #Load in metacyc rxn information, save those with associated KEGG ids
    metacyc_rxns = fread("metacyc/good_reactions.txt")
    metacyc_rxns = unique(metacyc_rxns[,list(`UNIQUE-ID`, Category, Value)])
    klinks = metacyc_rxns[Category=="DBLINKS" & grepl("LIGAND-RXN", Value)]
    klinks[,RID:=regmatches(Value, regexpr("R[0-9]+", Value))]
    klinks = klinks[RID %in% rxn_ids]
    metacyc_rxns = merge(metacyc_rxns, klinks[,list(`UNIQUE-ID`, RID)], by="UNIQUE-ID")
    setkey(metacyc_rxns, NULL)
    metacyc_rxns = metacyc_rxns[!duplicated(metacyc_rxns[,list(Category, Value, RID)])]
    #only save rxns with assoc metacyc info
    rxn_info = rxn_info[rxn_ids %in% metacyc_rxns[,RID]]
    rxn_ids = rxn_ids[rxn_ids %in% metacyc_rxns[,RID]] #this is a surprising drop
    rxn_kos = unique(unlist(sapply(rxn_info, function(x){ return( names(x$ORTHOLOGY))})))
    good_kos = kos[kos %in% rxn_kos]

    #Get compound ID info
    compound_key = fread("metacyc/good_compounds.txt")
    compound_key = unique(compound_key[,`UNIQUE-ID`, Value])
    setkey(compound_key, `UNIQUE-ID`)
    compound_key = compound_key[`UNIQUE-ID` %in% metacyc_rxns[,Value]]
    klinks = compound_key[grepl("LIGAND-CPD",Value)]
    #Get linked KEGG Compound IDs
    klinks[,CID:=toupper(substr(Value, 15, 20))]
    compound_key = unique(klinks[,list(`UNIQUE-ID`, CID)])
    setnames(compound_key, c("Value", "CID"))
    compounds = sort(unique(compound_key[,CID]))

    ##get stoichiometry and directions from metacyc
    stoich_mat = matrix(rep(NA), nrow = length(compounds), ncol = length(good_kos))
    net_rxns = c()
    net_kos = c()
    net_prods = c()
    net_reacs = c()
    net_dir = c()
    for(j in 1:length(rxn_ids)){ #get detailed info on each rxn
      rxn = rxn_ids[j]
      meta_info = metacyc_rxns[RID==rxn_ids[j]]
      meta_compounds = meta_info[Category %in% c("LEFT","RIGHT")]
      meta_compounds = merge(meta_compounds, compound_key, all.x=T, all.y=F, by="Value")[!is.na(CID)]
      #get reactants and products
      kos_involved = names(rxn_info[[j]]$ORTHOLOGY)
      kos_involved = kos_involved[kos_involved %in% good_kos]
      ko_inds = match(kos_involved,good_kos)
      reac = meta_compounds[Category=="LEFT", CID]
      prod = meta_compounds[Category=="RIGHT", CID]
      reac = reac[(!(reac %in% prod)) & (reac %in% compounds)] #get rid of things on both sides of equation and only save interesting compounds
      prod = prod[(!(prod %in% reac)) & (prod %in% compounds)]
      if(length(reac)==0 | length(prod)==0) next
      coefs = rep(1, length(reac)) #get stoichiometry
      for(k in 1:length(reac)){
        meta_id = meta_compounds[CID==reac[k], Value]
        ind = which(meta_info[,Value]==meta_id & meta_info[,Category=="LEFT"]) #coefficient always follows in dataset
        if(meta_info[ind+1,Category]=="^COEFFICIENT") coefs[k] = as.numeric(meta_info[ind+1, Value])
      }
      coefs2 = rep(1, length(prod))
      for(k in 1:length(prod)){
        meta_id = meta_compounds[CID==prod[k], Value]
        ind = which(meta_info[,Value]==meta_id & meta_info[,Category]=="RIGHT")
        if(ind < dim(meta_info)[1]){
          if(meta_info[ind+1,Category]=="^COEFFICIENT") coefs2[k] = as.numeric(meta_info[ind+1, Value])
        }
      }
      #get direction
      direction = meta_info[Category=="REACTION-DIRECTION", Value]
      if(length(direction)==0) dir = 0 else if(all(grepl("LEFT-TO-RIGHT",direction))) dir = 1 else if(direction=="REVERSIBLE") dir = 0 else if(all(grepl("RIGHT-TO-LEFT",direction))) dir = -1
      #still saving reversible rxn links
      for(k in 1:length(reac)){
        cmpd = grep(reac[k],compounds)
        if(length(cmpd) > 0){
          #slots that were NA before
          stoich_mat[cmpd, ko_inds][is.na(stoich_mat[cmpd, ko_inds])] = -1*dir*coefs[k]
          #non-NA slots
          stoich_mat[cmpd, ko_inds][!is.na(stoich_mat[cmpd, ko_inds])] = stoich_mat[cmpd, ko_inds][!is.na(stoich_mat[cmpd, ko_inds])] - 1*dir*coefs[k]
        }
      }
      #repeat for products
      for(k in 1:length(prod)){
        cmpd = grep(prod[k],compounds)
        if(length(cmpd) > 0){
          #slots that were NA before
          stoich_mat[cmpd, ko_inds][is.na(stoich_mat[cmpd, ko_inds])] = dir*coefs2[k]
          #non-NA slots
          stoich_mat[cmpd, ko_inds][!is.na(stoich_mat[cmpd, ko_inds])] = stoich_mat[cmpd, ko_inds][!is.na(stoich_mat[cmpd, ko_inds])] + dir*coefs2[k]
        }
      }
      n_subrxns = length(reac)*length(prod)*length(kos_involved)
      net_rxns = c(net_rxns, rep(rxn, n_subrxns))
      net_dir = c(net_dir, rep(dir, n_subrxns))
      for(i in 1:length(kos_involved)){
        for(p in 1:length(reac)){
          for(m in 1:length(prod)){
            net_kos = c(net_kos, kos_involved[i])
            net_reacs = c(net_reacs, reac[p])
            net_prods = c(net_prods, prod[m])
          }
        }
      }
    }
    network_table = data.table(Rxn = net_rxns, KO = net_kos, Reac = net_reacs, Prod = net_prods, Direction = net_dir)
    network_table = network_table[Prod != Reac]
    degree = apply(stoich_mat, 1, function(x){ length(x[!is.na(x)])})
    compounds = compounds[degree != 0]
    stoich_mat = stoich_mat[degree != 0,]
    degree = degree[degree != 0]
    if(degree_filter != 0){
      cat("Filtering currency metabolites\n")
      stoich_mat = stoich_mat[degree < degree_filter,]
      compounds = compounds[degree < degree_filter]
      degree = degree[degree < degree_filter]
    }
    degree = apply(stoich_mat, 1, function(x){ length(x[!is.na(x)])}) #double check
    compounds = compounds[degree != 0]
    stoich_mat = stoich_mat[degree != 0,]
    network_table = network_table[Reac %in% compounds & Prod %in% compounds]

    network_mat = stoich_mat
    network_mat[is.na(network_mat)] = 0
    if(normalize){
      negsums = apply(network_mat, 1, function(x){ sum(x[x < 0])})
      possums = apply(network_mat, 1, function(x){ sum(x[x > 0])})
      metlen = dim(network_mat)[1]
      for(j in 1:metlen){
        network_mat[j,][network_mat[j,] < 0] = network_mat[j,][network_mat[j,] < 0]/abs(negsums[j])
        network_mat[j,][network_mat[j,] > 0] = network_mat[j,][network_mat[j,] > 0]/possums[j]
      }
    }
    network_mat = data.frame(network_mat)
    names(network_mat) = good_kos
    row.names(network_mat) = compounds
    return(list(network_mat, stoich_mat, network_table))
  }
}

#' Calculate CMP scores based on community network and gene abundances
#'
#' @import data.table
#' @param emm Stoichiometric network matrix produced by generate_genomic_network
#' @param norm_kos Data.table of gene abundances
#' @return Data.table of CMP scores
#' @examples
#' get_prmt_scores(ko_net[[1]], gene_abunds)
#' @export
get_prmt_scores = function(emm, norm_kos){
  metlen = dim(emm)[1]
  nsamp = dim(norm_kos)[2]-1
  norm_kos_sub = norm_kos[KO %in% names(emm)]
  subjects = names(norm_kos)[names(norm_kos)!="KO"]
  prmt = matrix(rep(NA),nrow = metlen,ncol = nsamp)
  emm = emm[,match(norm_kos_sub[,KO], names(emm))]
  if(all(names(emm)==norm_kos_sub[,KO])){
    for(m in 1:nsamp){
      prmt[,m] =  as.matrix(emm) %*% unlist(norm_kos_sub[,subjects[m],with=F])
    }
  }else if(all(sort(names(emm))==sort(norm_kos_sub[,KO]))){ #Just out of order
    emm = emm[,order(names(emm))]
    norm_kos_sub = norm_kos_sub[order(KO)]
  } else {
    stop("Double check you are using the correct metabolic network! Gene KOs do not equal network KOs")
  }
  prmt = data.table(prmt,row.names(emm))
  setnames(prmt,c(subjects,"compound"))
  setkey(prmt,compound)
  return(prmt)
}

#' Make a data vector into a pairwise difference matrix (can be used for either PRMT scores or metabolite concentrations)
#'
#' @import data.table
#' @param metabolite metabolite ID
#' @param met_mat Abundances or scores for that metabolite across samples
#' @param diff_function Difference by default, could be fold_change
#' @return Data.frame of score or abundance pairwise differentials
#' @examples
#'
#' @export
make_pairwise_met_matrix = function(metabolite, met_mat, diff_function = "difference"){
  #diff_function options: 'difference', 'fold_change'
  #met_mat can be either PRMT scores or metabolomic abundances
  nsamp = ncol(met_mat) - 1
  subjects = names(met_mat)[1:nsamp]
  met_mat_small = unlist(met_mat[metabolite,subjects,with=F])
  met_matrix = matrix(rep(NA),nrow = nsamp,ncol = nsamp)
  if(diff_function == "difference"){
    for(j in 1:nsamp){
      for(k in 1:nsamp){
        met_matrix[j,k] = met_mat_small[k] - met_mat_small[j]
      }
    }
  }
  met_matrix = data.frame(met_matrix)
  names(met_matrix) = subjects
  row.names(met_matrix) = subjects
  return(met_matrix)
}

#' Modification of vegan mantel test to test 2-sided or significantly less than
#'
#' @import data.table
#' @param xdis A distance matrix
#' @param ydis Another distance matrix
#' @param method Correlation coefficient to use (pearson or spearman)
#' @param permutations Number of permutations to perform
#' @return Data.frame of score or abundance pairwise differentials
#' @examples
#'
#' @export
mantel_2sided = function (xdis, ydis, method = "pearson", permutations = 999,
                          strata, na.rm = FALSE, direction = "pos") {
  #direction must be either "pos","neg", or "two-sided"
  xdis <- as.dist(xdis)
  ydis <- as.vector(as.dist(ydis))
  if (na.rm)
    use <- "complete.obs"
  else use <- "all.obs"
  statistic <- cor(as.vector(xdis), ydis, method = method,
                   use = use)
  variant <- match.arg(method, eval(formals(cor)$method))
  variant <- switch(variant, pearson = "Pearson's product-moment correlation",
                    kendall = "Kendall's rank correlation tau", spearman = "Spearman's rank correlation rho",
                    variant)
  N <- attr(xdis, "Size")
  if (length(permutations) == 1) {
    if (permutations > 0) {
      arg <- if (missing(strata))
        NULL
      else strata
      permat <- t(replicate(permutations, permute::shuffle(N)))
      #permat <- t(replicate(permutations, permuted.index(N, strata = arg)))
    }
  }
  else {
    permat <- as.matrix(permutations)
    if (ncol(permat) != N)
      stop(gettextf("'permutations' have %d columns, but data have %d observations",
                    ncol(permat), N))
    permutations <- nrow(permutations)
  }
  if (permutations) {
    perm <- numeric(permutations)
    xmat <- as.matrix(xdis)
    asdist <- row(xmat) > col(xmat)
    ptest <- function(take, ...) {
      permvec <- (xmat[take, take])[asdist]
      drop(cor(permvec, ydis, method = method, use = use))
    }
    perm <- sapply(1:permutations, function(i, ...) ptest(permat[i, ], ...))
    if(direction=="pos"){
      signif <- (sum(perm >= statistic) + 1)/(permutations + 1)
    } else if(direction=="neg"){
      signif = (sum(perm <= statistic) + 1)/(permutations + 1)
    } else{
      signif = (sum(perm >= statistic) + sum(perm <= statistic) +1)/(permutations +1)
    }
  }
  else {
    signif <- NA
    perm <- NULL
  }
  res <- list(call = match.call(), method = variant, statistic = statistic,
              signif = signif, perm = perm, permutations = permutations)
  if (!missing(strata)) {
    res$strata <- deparse(substitute(strata))
    res$stratum.values <- strata
  }
  class(res) <- "mantel"
  return(res)
}

select_best_id2 = function(met_table, met_data, net_compounds, final_method = "first"){ ###no retention time, swedish data format
  met_table2 = data.frame(met_table)
  met_table2$Mass = as.character(met_table2$Mass)
  met_data = data.frame(met_data)
  subjects = names(met_data)[names(met_data)!="Mass"]
  met_data$Mass = as.character(met_data$Mass)
  met_table2 = met_table2[met_table2$KEGG %in% net_compounds,]
  met_data = met_data[met_data$Mass %in% met_table2$Mass,]
  met_table2 = met_table2[met_table2$Mass %in% met_data$Mass,]
  met_table2 = unique(met_table2[order(met_table2$Mass),])
  unique_ids = unique(unlist(met_table2$KEGG))
  unique_ids = sort(unique_ids[!is.na(unique_ids)]) #get list of potential IDs
  new_mets = head(met_data[,names(met_data)!="Mass"],0)
  good_mets = c()
  for(j in 1:length(unique_ids)){
    foo = met_table2[apply(met_table2,1,function(x){ unique_ids[j] %in% x}),]
    foo_good = foo[which.min(foo$Delta),]
    foo_data = met_data[met_data$Mass %in% foo$Mass,]
    #foo_good_close = foo[abs(foo$RetTime - foo_good$RetTime)<0.1,]
    if(dim(foo_data)[1] == 1){
      new_mets = rbind(new_mets, foo_data[,subjects])
      #figure out way to merge ones that are close in mass
      good_mets = c(good_mets, unique_ids[j])
      met_data = met_data[-which(met_data$Mass %in% foo$Mass),] #remove these peaks from consideration for future identifications
      met_table2 = met_table2[-which(met_table2$Mass %in% foo$Mass),]
    }
    if(dim(foo_data)[1] > 1){
      test = apply(foo_data[,1:12],2,function(x){length(x[x!=0])})
      if(length(test[test>1])>2){
        foo_good=foo[which.min(foo$Delta),]
        foo_data2=foo_data[foo_data$Mass %in% foo_good$Mass,subjects]
      } else foo_data2=apply(foo_data[,1:12],2,sum)
      new_mets = rbind(new_mets, foo_data2)
      good_mets = c(good_mets,unique_ids[j])
      met_data = met_data[-which(met_data$Mass %in% foo$Mass),] #remove these peaks from consideration for future identifications
      met_table2 = met_table2[-which(met_table2$Mass %in% foo$Mass),]
    }
  }
  #for(k in 1:length(unique_ids))
  #  new_mets = rbind(new_mets, lapply(met_data[which(single_final == unique_ids[k]),],sum))
  new_mets = data.table(new_mets, KEGG = good_mets)
  setkey(new_mets, KEGG)
  return(new_mets)
}

#select_best_id(kegg_cecum, cecum_mets_mass_numeric, net_compounds)
select_best_id = function(met_table, met_data, net_compounds, final_method = "first"){
  #takes list of compounds from a network, preferentially selects IDs in that network from putative options,
  #and combines data assigned to the same id
  #final_method can be 'first', 'random', or 'most_rxns' - how to make the final selection between multiple options
  #most_rxns not yet implemented but i think worth doing - forgot about this
  met_table2 = met_table[met_table$Delta!="-",]
  met_data = data.frame(met_data)
#  met_table3 = met_table2[apply(met_table2[,3:10],1,function(x){ length(x[!is.na(x)])})==1,]
  #only save IDs out of network if there are none in the network for that ion
   for(x in 1:length(met_table2[,1])){
     if(any(met_table2[x,3:length(met_table2[1,])] %in% net_compounds))
       met_table2[x,3:length(met_table2[1,])][!(met_table2[x,3:length(met_table2[1,])] %in% net_compounds)] = NA
   }
  met_table2 = met_table2[apply(met_table2,1,function(x){ any(!is.na(x[3:length(x)]))}),]
  for(j in 1:dim(met_table2)[1]){
    foo = unlist(sort(met_table2[j,3:dim(met_table2)[2]]))
    met_table2[j,3:dim(met_table2)[2]] = c(foo, rep(NA, dim(met_table2)[2]-2-length(foo)))
  }
  met_table2 = met_table2[order(met_table2$V1),]
  unique_ids = unique(unlist(met_table2[3:dim(met_table2)[2]]))
  unique_ids = sort(unique_ids[!is.na(unique_ids)])
#   dims = sapply(1:length(unique_ids), function(x){
#     foo= met_table2[apply(met_table2,1,function(y){ any(y[3:length(y)]==unique_ids[x] & !is.na(y[3:length(y)])) }),]
#     return(dim(foo)[1])
#   })
  new_mets = head(met_data[,names(met_data)!="Mass"],0)
  good_mets = c()
  for(j in 1:length(unique_ids)){
    foo = met_table2[apply(met_table2,1,function(x){ unique_ids[j] %in% x}),]
    foo_good = foo[which.min(foo$Delta),]
    foo_data = met_data[c(as.numeric(row.names(foo))),]
    foo_good_close = foo_good[abs(foo$RetTime - foo_good$RetTime)<0.1,]
    if(dim(foo_good_close)[1] > 0){
      new_mets = rbind(new_mets, lapply(met_data[row.names(met_data) %in% as.numeric(row.names(foo_good_close)),names(new_mets)],sum))
      good_mets = c(good_mets, unique_ids[j])
      met_data = met_data[!(row.names(met_data) %in% row.names(foo_good_close)),] #remove these peaks from consideration for future identifications
      met_table2 = met_table2[!(row.names(met_table2) %in% row.names(foo_good_close)),]
    }
  }
  #add together matching ions
  #for(k in 1:length(unique_ids))
  #  new_mets = rbind(new_mets, lapply(met_data[which(single_final == unique_ids[k]),],sum))
  new_mets = data.table(new_mets, KEGG = good_mets)
  setkey(new_mets, KEGG)
  return(new_mets)
}

#' Runall function to do complete predictions and comparison for all shared metabolites
#'
#' @import data.table
#' @param genes Gene abundances
#' @param mets Metabolite abundances
#' @param file_prefix Prefix for output files
#' @param correction Type of multiple hypothesis correction to perform (bonferroni or fdr)
#' @param cutoff Q/P-value cutoff for significance
#' @param net_method Network generation method (see generate_genomic_network)
#' @param id_met Whether metabolites have putative identifications that need to be processed
#' @param met_id_file If id_met, file of metabolite identifications
#' @param degree_filter Threshold for filtering currency metabolites
#' @param minpath_file
#' @param cor_method
#' @param net_file
#' @param nperm
#' @param nonzero_filter
#' @return No return, writes output to file
#' @examples
#'
#' @export
run_all_metabolites = function(genes, mets, file_prefix = 'net1', correction = "fdr", cutoff = 0.1,
                               net_method = "load", id_met = F, met_id_file = '',
                               degree_filter = 0, minpath_file = '', cor_method = "spearman",
                               net_file = "", nperm = 20000, nonzero_filter=2){
  #must be data.tables keyed by KOs/KEGG IDs in first column, all other columns are subject IDs
  #correction must be either "bonferroni" or "fdr", cutoff is q value cutoff
  #id_mets specifies whether to use network for improved metabolite identification (i.e. for Braun datasets)
  #if true, text file with list of putative met ids from Metabosearch must be included
  kos = genes[,KO]
  nsamp = dim(genes)[2] - 1
  subjects = intersect(names(genes), names(mets))
  ko_net = generate_genomic_network(kos, keggSource = net_method, degree_filter = degree_filter, minpath_file = minpath_file, networkFile = net_file)
  ko_network_mat = ko_net[[1]]
  stoich_mat = ko_net[[2]]
  ko_net_table = ko_net[[3]]
  if(id_met){
    cat("Selecting best metabolite identifications\n")
    net_compounds = row.names(ko_net[[1]])
    met_id_table = read.table(met_id_file, sep = "\t", header=T)
    met_id_table = met_id_table[met_id_table$Delta!="-",] #get rid of mets with no possible KEGG ID
    #met_id_list = strsplit(unlist(met_id_list), split = " ")
    if(grepl("swedish",file_prefix)){
      new_mets = select_best_id2(met_id_table, mets, net_compounds)
      cat("Swedish ID met\n")
    } else new_mets = select_best_id(met_id_table, mets, net_compounds)
    mets = new_mets
    cat("Done identifying metabolites!\n")
  }

  #get mets
  metIDs = mets[,KEGG]
  shared_mets = metIDs[metIDs %in% row.names(ko_network_mat)]
  #get PRMT scores
  norm_kos = genes[,c(subjects,"KO"),with=F]
  prmt_mat = get_prmt_scores(ko_network_mat, norm_kos)

  #do all comparisons
  all_comparisons = vector("list",length(shared_mets))
  for(j in 1:length(shared_mets)){
    good_subs = intersect(names(mets)[which(!is.na(unlist(mets[shared_mets[j],subjects,with=F])))], names(prmt_mat)[which(!is.na(unlist(prmt_mat[shared_mets[j],subjects,with=F])))])
    prmt_vector = unlist(prmt_mat[shared_mets[j],good_subs,with=F])
    met_vector = unlist(mets[shared_mets[j],good_subs,with=F])
    #check for too many 0s
    if(length(met_vector[met_vector!=0]) <= nonzero_filter | length(prmt_vector[prmt_vector!=0]) <= nonzero_filter){
      all_comparisons[[j]] = NA
    }else{
      met_mat = make_pairwise_met_matrix(shared_mets[j], prmt_mat[,c(good_subs, "compound"),with=F])
      metabol_mat = make_pairwise_met_matrix(shared_mets[j], mets[,c(good_subs,"KEGG"),with=F])
      test = vegan::mantel(met_mat,metabol_mat,method=cor_method,permutations = nperm)
      test_n = mantel_2sided(met_mat,metabol_mat,method=cor_method,permutations = nperm,
                             direction = "neg")
      all_comparisons[[j]] = list(ID = shared_mets[j], PRMT = met_mat, Mets = metabol_mat, Mantel = list(test,test_n))
    }
  }

  #remove ones that didn't pass filter
  shared_mets = shared_mets[which(!is.na(all_comparisons))]
  all_comparisons = all_comparisons[which(!is.na(all_comparisons))]

  #whole set analysis
  cors_s = sapply(all_comparisons,function(x){return(x$Mantel[[1]]$statistic)})
  pvals_s = sapply(all_comparisons,function(x){return(x$Mantel[[1]]$signif)})
  pvals2_s = correct(pvals_s, method = correction)

  cors_n = sapply(all_comparisons,function(x){return(x$Mantel[[2]]$statistic)})
  pvals_n = sapply(all_comparisons,function(x){return(x$Mantel[[2]]$signif)})
  pvals2_n = correct(pvals_n, method = correction)

  node_data = data.table(compound = shared_mets, CorrS = cors_s, PValS = pvals_s, QValS = pvals2_s,
                         CorrN = cors_n, PValN = pvals_n, QValN = pvals2_n)
  setkey(node_data,compound)

  #save everything
  #write edge file
  write.table(ko_net_table,file=paste(file_prefix,'_edges.txt',sep=''),sep="\t",quote=F,row.names=F)
  #write node attribute file
  write.table(node_data,file = paste(file_prefix,'_nodes.txt',sep=''),sep="\t",quote=F,row.names=F)

  signif_pos = node_data[!is.na(node_data$QValS) & (node_data$QValS < cutoff) & node_data$PValS < 0.05,]
  signif_neg = node_data[!is.na(node_data$QValN) & node_data$QValN < cutoff & node_data$PValN < 0.05,]
  write.table(signif_pos, file = paste(file_prefix,'_signifPos.txt',sep = ''), sep = "\t", quote = F, row.names = F)
  write.table(signif_neg, file = paste(file_prefix,'_signifNeg.txt',sep = ''), sep = "\t", quote = F, row.names = F)
  save(norm_kos, mets, ko_net, all_comparisons, node_data, file = paste(file_prefix,'_out.rda',sep=''))
  return(NULL)
}


#' Runall function to do complete predictions and evaluate classification of metabolites as high or low abundance
#'
#' @import data.table
#' @param genes Gene abundances
#' @param mets Metabolite abundances
#' @param file_prefix Prefix for output files
#' @param correction Type of multiple hypothesis correction to perform (bonferroni or fdr)
#' @param cutoff Q/P-value cutoff for significance
#' @param net_method Network generation method (see generate_genomic_network)
#' @param id_met Whether metabolites have putative identifications that need to be processed
#' @param met_id_file If id_met, file of metabolite identifications
#' @param degree_filter Threshold for filtering currency metabolites
#' @param minpath_file
#' @param net_file
#' @param quant
#' @param plot_rank
#' @param plot_continuous
#' @param nonzero_filter
#' @return No return, writes output to file
#' @examples
#'
#' @export
run_all_metabolites2 = function(genes, mets, file_prefix = 'net1', correction = "fdr", cutoff = 0.1,
                               net_method = "load", id_met = F, met_id_file = '',
                               degree_filter = 0, minpath_file = '',
                               net_file = "", quant = 0.5, plot_rank = F, plot_continuous=T, nonzero_filter=2, rel_abund_mets = F){
  #must be data.tables keyed by KOs/KEGG IDs in first column, all other columns are subject IDs
  #correction must be either "bonferroni" or "fdr", cutoff is q value cutoff
  #id_mets specifies whether to use network for improved metabolite identification (i.e. for Braun datasets)
  #if true, text file with list of putative met ids from Metabosearch must be included
  #can predict above/below median or other threshold with quant argument
  kos = genes[,KO]
  nsamp = dim(genes)[2] - 1
  subjects = intersect(names(genes), names(mets))
  ko_net = generate_genomic_network(kos, keggSource = net_method, degree_filter = degree_filter, minpath_file = minpath_file,
                                    networkFile = net_file)
  ko_network_mat = ko_net[[1]]
  stoich_mat = ko_net[[2]]
  ko_net_table = ko_net[[3]]
  if(id_met){
    cat("Selecting best metabolite identifications\n")
    net_compounds = row.names(ko_net[[1]])
    met_id_table = read.table(met_id_file, sep = "\t", header=T)
    #met_id_table = met_id_table[met_id_table$Delta!="-",] #get rid of mets with no possible KEGG ID
    #met_id_list = strsplit(unlist(met_id_list), split = " ")
    if(grepl("swedish",file_prefix)){
      new_mets = select_best_id2(met_id_table, mets, net_compounds)
    } else new_mets = select_best_id(met_id_table, mets, net_compounds)
    mets = new_mets
    cat("Done identifying metabolites!")
  }
  metIDs = mets[,KEGG]
  norm_kos = normalize_ko_counts(genes,subjects, logdata = F) #this actually does nothing right now
  #norm_kos = norm_kos[KO %in% names(ko_network_mat)]
  prmt_mat = get_prmt_scores(ko_network_mat, norm_kos)
  shared_mets = metIDs[metIDs %in% row.names(ko_network_mat)]
  all_comparisons = vector("list",length(shared_mets))
  #The idea here is to simply classify samples as higher or lower than the median and just look at classification error
  #rather than Mantel test
  for(j in 1:length(shared_mets)){
    prmts = unlist(prmt_mat[shared_mets[j],subjects,with=F])
    met1 = unlist(mets[shared_mets[j],subjects, with=F]) #these have very different distributions - I wonder how often this is the case
    if(length(met1[met1!=0]) <= nonzero_filter | length(prmts[prmts!=0]) <= nonzero_filter){
      all_comparisons[[j]] = NA
    }else{
      med_met = quantile(met1, probs = quant)
      if(any(met1==med_met)) med_prmt = median(prmts[met1==med_met]) else med_prmt = mean(prmts[which(abs(met1-med_met) <= min(abs(met1-med_met)) + 0.001)]) #median of two closest from which median was calculated
      higher_pred = ifelse(prmts > med_prmt,1,0)
      higher_obs = ifelse(met1 > med_met, 1, 0)
      error = sum(abs(higher_pred-higher_obs))/length(higher_obs)
      #need to clarify - this is different from not having any information
      if(plot_rank){
        requireNamespace("ggplot2", quietly = TRUE)
        preds = data.frame(Prediction = factor(higher_pred), Prmt = prmts, Value = met1, Rank = rank(met1))
        ggplot2::ggplot(preds, ggplot2::aes(x=Rank,y=Value, col = Prediction)) + ggplot2::geom_point(size=3) + ggplot2::xlim(c(0,max(preds$Rank)))+ ggplot2::theme_bw() + ggplot2::scale_color_manual(values=c("purple","orange")) +
          ggplot2::annotate("text",x=0.3*max(preds$Rank), y=0.8*max(preds$Value),label=met_names(shared_mets[j]))+
          ggplot2::annotate("text",x=0.3*max(preds$Rank), y=0.6*max(preds$Value), label=paste(length(preds$Value[preds$Value==min(preds$Value) & preds$Prediction==0]), length(preds$Value[preds$Value==min(preds$Value) & preds$Prediction==1]), sep = "  "))
        ggplot2::ggsave(file=paste0(file_prefix,"_",shared_mets[j],".png"))
      }
      if(plot_continuous){
        if(any(prmts!=0)){
          plot_ref_mets_by_prmts(met1, prmts,shared_mets[j], file_prefix)
        }
      }
      if(any(prmts!=0)) {
        if(any(higher_pred != 0) & any(higher_pred != 1) & any(higher_obs !=0) & any(higher_obs !=1)){
          ftest = fisher.test(higher_pred, higher_obs, alternative = "greater")
          ftest2 = fisher.test(higher_pred, higher_obs, alternative = "less")
        } else {
          ftest = 1
          ftest2 = 1
        }
      }else {
        ftest = NA
        ftest2 = NA
      }
      all_comparisons[[j]] = list(ID = shared_mets[j], Error = error, Test = ftest, Test2 = ftest2, PRMT_med = med_prmt, Met_med = med_met, Pred=higher_pred, Obs=higher_obs)
    }
  }

  #remove ones that didn't pass filter
  shared_mets = shared_mets[which(!is.na(all_comparisons))]
  all_comparisons = all_comparisons[which(!is.na(all_comparisons))]

  #whole set analysis
  pvals = sapply(all_comparisons, function(x){ if(is.na(x$Test)) return(NA) else if(is.numeric(x$Test)) return(1) else return(x$Test$p.value)})
  pvals_c = correct(pvals, method = correction)
  sensitivity = sapply(all_comparisons, function(x){ length(x$Pred[x$Pred==x$Obs & x$Pred==1])/length(x$Obs[x$Obs==1])})
  specificity = sapply(all_comparisons, function(x){ length(x$Pred[x$Pred==x$Obs & x$Pred==0])/length(x$Obs[x$Obs==0])})
  precision = sapply(all_comparisons, function(x){ length(x$Pred[x$Pred==1 & x$Pred==x$Obs])/length(x$Pred[x$Pred==1])})
  #cors_s = sapply(all_comparisons,function(x){return(x$Mantel[[1]]$statistic)})
  #pvals_s = sapply(all_comparisons,function(x){return(x$Mantel[[1]]$signif)})
  #pvals2_s = correct(pvals_s, method = correction)
  pvals_n = sapply(all_comparisons, function(x){if(is.na(x$Test)) return(NA) else if(is.numeric(x$Test2)) return(1) else return(x$Test2$p.value) })
  pvals_nc = correct(pvals_n, method = correction)
  accuracy = 1 - sapply(all_comparisons, function(x){ return(x$Error)})

  node_data = data.table(compound = shared_mets, PValS = pvals, QValS = pvals_c,
                         #CorrP = cors_p, PValP = pvals_p, QValP = pvals2_p,
                         PValN = pvals_n, QValN = pvals_nc, Accuracy = accuracy, Sensitivity = sensitivity, Specificity = specificity, Precision = precision)
  setkey(node_data,compound)
  #write to network file
  write.table(ko_net_table,file=paste(file_prefix,'_edges.txt',sep=''),sep="\t",quote=F,row.names=F)
  #write node attribute file
  write.table(node_data,file = paste(file_prefix,'_nodes2.txt',sep=''),sep="\t",quote=F,row.names=F)
  #plot_mantel_results(all_comparisons, node_data, file_prefix)
  signif_pos = node_data[!is.na(node_data$QValS) & (node_data$QValS < cutoff),]
  signif_neg = node_data[!is.na(node_data$QValN) & node_data$QValN < cutoff,]
  #write.table(signif_pos, file = paste(file_prefix,'_signifPos.txt',sep = ''), sep = "\t", quote = F, row.names = F)
  #write.table(signif_neg, file = paste(file_prefix,'_signifNeg.txt',sep = ''), sep = "\t", quote = F, row.names = F)
  save(norm_kos, mets, ko_net, all_comparisons, node_data, file = paste(file_prefix,'_out2.rda',sep=''))
  return(NULL)
  #return(list(all_comparisons = all_comparisons, met_table = node_data))
}

#' Multiple hypothesis correction
#'
#' @import data.table
#' @param pvals Vector of p-values
#' @param method Must be either "fdr" or "bonferroni"
#' @return Vector of corrected values
#' @examples
#'
#' @export
correct = function(pvals, method = "fdr"){
 #vector of p values, method must be either "bonferroni" or "fdr"
  foo = which(!is.na(pvals))
  if(method == "bonferroni") {
    pvalsAdj = p.adjust(pvals[!is.na(pvals)], method="bonferroni")
    pvals2 = sapply(1:length(pvals), function(x){ if(is.na(pvals[x])) return(NA) else return(pvalsAdj[which(foo==x)])})
  }else{
    qvals = qvalue::qvalue(pvals[foo])$qvalues
    pvals2 = sapply(1:length(pvals), function(x){ if(is.na(pvals[x])) return(NA) else return(qvals[which(foo==x)])})
  }
}

plot_mantel_results = function(metabolite_list, node_data, file_prefix){
  pdf(file = paste(file_prefix,"_qval_hist.pdf", sep = ''), width = 10)
  par(mfrow = c(1,3))
  hist(node_data$QValS, main = "Spearman correlation Mantel", ylim = c(0, 0.6*length(metabolite_list)), breaks = 25)
  abline(v = 0.05, col = "red")
  #hist(node_data$QValP, main = "Pearson correlation Mantel", ylim = c(0, 0.6*length(metabolite_list)), breaks = 25)
  #abline(v = 0.05, col = "red")
  hist(node_data$QValN, main = "Negative correlation Mantel", ylim = c(0, 0.6*length(metabolite_list)), breaks = 25)
  abline(v = 0.05, col = "red")
  dev.off()
}


#' Get metabolite name for KEGG ID(s)
#'
#' Requires metaboliteCategories_compoundNames.txt to be in the same directory
#'
#' @import data.table
#' @param met_ids One or a vector of KEGG compound IDs
#' @return vector of character metabolite names
#' @examples
#' met_names(c("C00082", "C00334"))
#'
#' @export
met_names = function(met_ids){
  #path_key = fread("metaboliteCategories_compoundNames.txt")
  setkey(path_key, "compound")
  return(path_key[met_ids, CompoundName])
#   load("KeggCompoundNames.rda")
#   return(sapply(met_ids, function(x){ return(all_met_names[match(x,names(all_met_names))][[1]][1])}))
}


plot_ref_mets_by_prmts = function(met, prmts, id, file_prefix){
  requireNamespace("ggplot2", quietly = TRUE)
  med_met = quantile(met, probs = quant)
  if(any(met==med_met)) med_prmt = median(prmts[met==med_met]) else med_prmt = mean(prmts[which(abs(met-med_met) <= min(abs(met-med_met)) + 0.001)]) #median of two closest from which median was calculated
  ref_met = met - med_met
  ref_prmt = prmts - med_prmt
  higher_pred = ifelse(prmts > med_prmt,1,0)
  higher_obs = ifelse(met > med_met, 1, 0)
  plot_data = data.frame(PRMT=ref_prmt, Met=ref_met, HigherPred=factor(higher_pred),HigherObs=factor(higher_obs))
  x0=ggplot2::ggplot(plot_data,ggplot2::aes(x=PRMT,y=Met))+ggplot2::geom_point(size=4)+ggplot2::xlab("Score reference difference")+ ggplot2::ylab("Metabolite reference difference") + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept=med_prmt, linetype=2) + ggplot2::geom_hline(ylintercept=med_met, linetype=2)#+
  print(med_met)
  print(med_prmt)
    #annotate("text",x=0.3*max(plot_data$PRMT), y=0.8*max(plot_data$Met),label=met_names(id))
  return(x0)
  #ggsave(file=paste0(file_prefix,id,".png"))
}


ggMMplot <- function(var1, var2, prop=T, fontsize=7, text_location="bottom"){
  requireNamespace("ggplot2", quietly = TRUE)
  levVar1 <- length(levels(var1))
  levVar2 <- length(levels(var2))
  y = ifelse(text_location=="bottom",-0.1,1.1)
  ylims = ifelse(text_location=="bottom", c(-0.1,1),c(0,1.1))
  if(prop) jointTable <- prop.table(table(var1, var2)) else jointTable = table(var1, var2)
  plotData <- as.data.frame(jointTable)
  if(prop) plotData$marginVar1 <- prop.table(table(var1)) else plotData$marginVar1 = table(var1)
  plotData$var2Height <- plotData$Freq / plotData$marginVar1
  plotData$var1Center <- c(0, cumsum(plotData$marginVar1)[1:levVar1 -1]) + plotData$marginVar1 / 2

  tab1 = data.frame(table(var1))
  plotData$labels = ifelse(tab1[match(plotData$var1,plotData$var1),2]==0, "", paste0(plotData$var1, "\n(n=",tab1[match(plotData$var1,plotData$var1),2],")"))
  ggplot2::ggplot(plotData, ggplot2::aes(var1Center, var2Height)) +
    ggplot2::geom_bar(stat = "identity", ggplot2::aes(width = marginVar1, fill = var2), col = "Black") +
    ggplot2::annotate("text", label = plotData$labels, x = plotData$var1Center,y=y, size=fontsize) + ggplot2::scale_x_continuous(expand=c(0,0)) + ylim(ylims[1], ylims[2])
}

#' Get non-reverible reactions of network.
#'
#' @import data.table
#' @param pvals Vector of p-values
#' @param method Must be either "fdr" or "bonferroni"
#' @return Vector of corrected values
#' @examples
#'
#' @export
get_non_rev_rxns = function(rxn_table, all_rxns=T){ #whether to return all reactions or only 1/2 of each reversible reaction since info is redundant
  if(dim(rxn_table)[1]>0){
    all_sorted= data.table(t(apply(rxn_table[,list(KO,Reac,Prod,stoichReac,stoichProd)],1,function(y){ sort(unlist(y))})))
    all_sorted[,Count:=.N,by=names(all_sorted)]
    all_sorted[,Reversible:=ifelse(Count==2,1,0)]
    rxn_table[,Reversible:=ifelse(all_sorted[,Count]==2,1,0)] #order is the same
    if("V1" %in% names(all_sorted)){
      all_sorted = unique(all_sorted)
      all_sorted = all_sorted[,list(V5, V4, V3, V2, V1, Reversible)]
      setnames(all_sorted, c("KO","Reac", "Prod", "stoichReac", "stoichProd", "Reversible"))
    }
    if(all_rxns) return(rxn_table) else return(all_sorted)
  } else return(NULL)
}

#' Calculate scores based on only a single gene.
#' 
#' @import data.table
#' 
#' 
#' 
#' @export
single_gene_cmp = function(compound, gene, norm_kos, ko_net){
  
  
}

#' Identify potential important gene contributors for each metabolite.
#'
#' @import data.table
#' @param j
#' @param prmts_sub_good
#' @param all_rxns
#' @param subjects
#' @param norm_kos
#' @param ko_net
#' @return list of contributors and correlations
#' @examples
#'
#' @export
gene_contributions = function(j, prmts_sub_good, all_rxns, subjects, norm_kos, ko_net){
  #index (from sapply usually), prmt matrix, list of reaction tables for each compound, subjects, gene matrix, full network
  if(!is.null(all_rxns[[j]])){
    compound = prmts_sub_good[j,compound]
    kos_involved = unique(all_rxns[[j]][Reversible==0,KO])
    ko_vals = norm_kos[KO %in% kos_involved]
    #plan - look at correlation between old prmt score and prmt score w/out each KO - least correlated is most impactful KO
    ko_prmt_cors = sapply(kos_involved, function(x){
      vals_without = ko_vals[KO != x]
      prmt_without = data.matrix(ko_net[[1]][compound,vals_without[,KO]])%*%data.matrix(vals_without[,subjects,with=F])
      return(cor(as.vector(unlist(prmts_sub_good[compound,subjects,with=F])),as.vector(prmt_without), method="spearman"))
    })
    ko_cors = data.table(KO=kos_involved, Cor=ko_prmt_cors)
    #save KOs that have a major effect
    ko_good = ko_cors[is.na(Cor)|(Cor < 0.5),KO]
    net_primary = all_rxns[[j]][KO %in% ko_good & Reversible==0]
    if(dim(net_primary)[1]==0){ #if none have a major effect, look at all reactions to get primary predicting force
      if(all(all_rxns[[j]][Reversible==0,Prod]==compound)){ primary_make = 1} else if(all(all_rxns[[j]][Reversible==0,Reac]==compound)){ primary_make = -1} else primary_make = 0
    } else if(all(net_primary[,Prod]==compound)) {
      primary_make = 1
      } else if(all(net_primary[,Reac]==compound)) {
        primary_make = -1 } else { primary_make = 0 }
    nprod = length(unique(all_rxns[[j]][Reversible==0 & compound==Prod, KO]))
    nreac = length(unique(all_rxns[[j]][Reversible==0 & compound==Reac, KO]))
    if(nrow(net_primary)==0){
      nkey_prod = length(unique(all_rxns[[j]][Reversible==0 & compound==Prod, KO]))
      nkey_reac = length(unique(all_rxns[[j]][Reversible==0 & compound==Reac, KO]))
    } else{
      nkey_prod = length(unique(net_primary[Reversible==0 & compound==Prod, KO]))
      nkey_reac = length(unique(net_primary[Reversible==0 & compound==Reac, KO]))
    }
    return(list(ko_cors, net_primary, primary_make, c(nKOReac = nreac, nKOProd = nprod, nkeyKOReac = nkey_reac, nkeyKOProd = nkey_prod)))
    #  return(ko_cors)
  } else return(NULL)
}

#' Compare a single set each of CMP scores and metabolite concentrations
#'
#' @import data.table
#' @param met_met metabolite to use from metabolite concentration data
#' @param met_prmt metabolite to use from CMP score data
#' @param met_all matrix of metabolite concentrations
#' @param prmt_all matrix of CMP scores
#' @param posneg Whether to test for positive or negative concentration
#' @param cor_method spearman or pearson
#' @param nperm number of permutations for Mantel test
#' @return p-value from Mantel test
#' @examples
#' compare_met("C00001", "C00002", mets, cmp_mat, "pos")
#' @export
compare_met = function(met_met, met_prmt, met_all, prmt_all, posneg="pos", cor_method="spearman", nperm=15000){
  good_subs = names(met_all)[which(!is.na(unlist(met_all[met_met])) & names(met_all)!="KEGG")]
  met_mat = make_pairwise_met_matrix(met_prmt, prmt_all[,c(good_subs,"compound"),with=F])
  metabol_mat = make_pairwise_met_matrix(met_met, met_all[,c(good_subs,"KEGG"),with=F])
  if(posneg=="pos") test = vegan::mantel(met_mat,metabol_mat,method=cor_method,permutations = nperm)
  else test = mantel_2sided(met_mat,metabol_mat,method=cor_method,permutations = nperm, direction = "neg")
  return(test$signif)
}

#' Function for getting ko abundances from BV qPCR data
#'
#' @import data.table
#' @param qpcr qPCR species abundance data across samples
#' @param ref_kos gene content of every species
#' @param scale_factor Scale abundances by a factor
#' @return gene abundances across samples
#' @examples
#' kos_from_species(bv_qpcr, genome_content)
#' @export
kos_from_species = function(qpcr, ref_kos, scale_factor = 1){
  #kos_by_sample=list(NULL) ##kos by ksample and species - have qpcr for 14 species but genomes for only 11 of them
  ref_names = unique(ref_kos[,species])
  spec_names = names(qpcr)[!names(qpcr) %in% c("ID", "Sample")]
  spec_names = spec_names[!is.na(spec_names)]
  kolist=data.table()
  for(j in 1:nrow(qpcr)){
    for(i in spec_names){
      #ind=grep(names(qpcr)[i],ref_names,ignore.case=T)
      if(any(grepl(i, ref_names)) & !is.na(unlist(qpcr[j,i,with=F]))){
        amt=as.numeric(qpcr[j,i,with=F]) #amount in that sample
        if(amt>0){
          kolist=rbind(kolist, data.table(KO=as.character(ref_kos[species==i,KO]),Abund=ref_kos[species==i,CopyNum]*amt/scale_factor, Species=i, Sample=qpcr[j,ID]))
          ##scale by 1000 to make numbers more manageable
        }
      } #else kos_by_sample[[j]]=NULL
    }
  }
  BV_kos = kolist[,sum(Abund),by=list(KO,Sample)]
  BV_kos = dcast.data.table(BV_kos,KO~Sample, value.var="V1")
  return(BV_kos)
}

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
#' @return
#' @examples
#' sapply(1:length(metabolites), prmt_species_contributions, cmp_scores, all_rxns, all_subjects, ko_abunds, ko_net, spec_abunds, ref_kos)
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
#' sapply(1:length(metabolites), prmt_species_contributions, cmp_scores, all_rxns, all_subjects, ko_abunds, ko_net, spec_abunds, ref_kos)
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

#' Test against list of known microbial metabolites
#'
#' @import data.table
#' @param node_data output of run_all_metabolites
#' @param met_list List of metabolites to test against
#' @return Null, just prints a bunch of p-values
#' @examples
#' test_met_enrichment(node_data, microbial_mets)
#' @export
test_met_enrichment = function(node_data, met_list){
  node_data[,MicrobeControl:=ifelse(compound %in% met_list[,KEGG],1,0)]
  if(any(node_data[,MicrobeControl])) {
    print(table(node_data[,MicrobeControl],node_data[,Pos]|node_data[,Neg]))
    print(fisher.test(node_data[,MicrobeControl],node_data[,Pos],alternative="greater"))
    print(fisher.test(node_data[,MicrobeControl],node_data[,Neg], alternative="greater"))
    print(fisher.test(node_data[,MicrobeControl], node_data[,Pos]|node_data[,Neg], alternative="greater"))
    print(wilcox.test(node_data[MicrobeControl==0,CorrS],node_data[MicrobeControl==1,CorrS], alternative="greater"))
    print(wilcox.test(node_data[MicrobeControl==0,abs(CorrS)],node_data[MicrobeControl==1,abs(CorrS)], alternative="greater"))
  }
}

#Capitalize 1st letter of words
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}


#' Get network distances between 2 compounds
#'
#' @import data.table
#' @param c1 first compound
#' @param c2 second compound
#' @param allnet edge list output of generate_genomic_network
#' @param max_dist distance at which the function will stop looking and assume the two compounds are not connected
#' @return integer representing the number of reaction steps between the two compounds, or max_dist if they are not connected
#' @examples
#' get_net_dist("C00001", "C00002", allnet)
#' @export
get_net_dist = function(c1, c2, allnet, max_dist = 20){
  if(c1==c2) return(0)
  match_rxns = allnet[Prod==c1 | Reac==c1]
  netdist = 1
  while(!(c2 %in% match_rxns[,Prod]) & !(c2 %in% match_rxns[,Reac]) & netdist < max_dist){
    netdist = netdist + 1
    all_comps = unique(c(match_rxns[,Prod], match_rxns[,Reac]))
    match_rxns = allnet[Prod %in% all_comps | Reac %in% all_comps]
  }
  return(netdist)
}

# calculate_net_dist=function(compound1,compound2, netmat, maxdist){
#   #re-doing to make sure there's an actual reaction connecting the 2 compounds - net_mat should be edge list (network[[3]])
#   if(identical(compound1,compound2)) return(0)
#   level = 1
#   comps1 = sort(unique(c(netmat[Prod==compound1, Reac], netmat[Reac==compound1, Prod])))
#   while(!(compound2 %in% comps1) & level < 12){
#     level = level+1
#     comps1 = sort(unique(c(netmat[Prod %in% comps1, Reac], netmat[Reac %in% comps1, Prod])))
#   }
#   return(level)
# }
#

#' Defaults for NULL values
#' @export
`%||%` <- function(a, b) if (is.null(a)) b else a
