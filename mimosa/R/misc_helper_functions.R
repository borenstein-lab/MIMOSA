##Miscellaneous helper/utility functions

#' Read in processed files and assign keys to datasets. Sets NA metabolite values to 0.
#'
#' @import data.table
#' @param genefile File where gene data is located.
#' @param metfile File where metabolite data is located.
#' @return A list in which the first item is a data.table of gene abundances and the second is a data.table of metabolite abundances.
#' @examples
#' read_files(gene_file, met_file)
#' @export
#'
#' @useDynLib mimosa
#' @importFrom Rcpp sourceCpp
read_files = function(genefile, metfile){
  genes = fread(genefile, header=T, sep="\t")
  setkey(genes,KO)
  mets = fread(metfile, header=T, sep="\t")
  #save only samples that have both kinds of data and put datasets in the same order of subjects/samples
  subjects = sort(intersect(names(genes), names(mets)))
  if(length(subjects) < length(names(genes))-1 | length(subjects) < length(names(mets))-1){
    cat("Only using sample IDs found for both genes and metabolites\n")
  }
  genes = genes[,c("KO", subjects), with=F]
  if("KEGG" %in% names(mets)) mets = mets[,c(subjects,"KEGG"), with=F] else mets = mets[,c("Mass", subjects), with=F]
  if("KEGG" %in% names(mets)){
    mets = data.table(mets[,lapply(.SD, as.numeric), .SDcols = subjects], KEGG = mets[,KEGG])
    setkey(mets, KEGG)
  } #setkey(mets,KEGG) #2 possibilities for metabolite file format
  
  #Set characters to NAs, NAs to 0
  for(j in names(mets)[!names(mets) %in% c("KEGG", "Mass")]){
    set(mets,which(is.na(mets[[j]])),j,0)
  }
  return(list(genes, mets))
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

#' Get non-reverible reactions of network.
#'
#' @import data.table
#' @param pvals Vector of p-values
#' @param method Must be either "fdr" or "bonferroni"
#' @return Vector of corrected values
#' @examples
#' get_non_rev_rxns(rxn_table)
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


#' Convert edge list network format to matrix format, including row normalization.
#'
#' @import data.table
#' @param rxn_table Edge list format of network, as from 3rd output of generate_genomic_network
#' @return list of two matrices, one with NAs filled in with zeros, one without, as first two outputs of generate_genomic_network
#' @examples
#' make_network_matrix(randomized_edge_list)
#' @export
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
  if(all(as.vector(xdis) == 0)|all(ydis == 0)){
    stop("All values are 0 for one distance matrix")
  }
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
      permat = make_perm_mat(N, permutations)
      #permat <- t(replicate(permutations, permute::shuffle(N)))
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

#' Get metabolite name for KEGG ID(s)
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

#' Function for getting ko abundances from species abundance data and table of reference genomes
#'
#' @import data.table
#' @param spec_abunds Species abundance data table across samples
#' @param ref_kos gene content of every species
#' @param scale_factor Scale abundances by a factor
#' @return gene abundances across samples
#' @examples
#' kos_from_species(bv_qpcr, genome_content)
#' @export
kos_from_species = function(spec_abunds, ref_kos, scale_factor = 1){
  ref_names = unique(ref_kos[,Species])
  spec_names = spec_abunds[,unique(Species)]
  spec_names = spec_names[!is.na(spec_names)]
  samp_names = spec_abunds[,unique(Sample)]
  kolist=data.table()
  for(j in 1:length(samp_names)){
    for(i in spec_names){
      if(i %in% ref_names & nrow(spec_abunds[Sample==samp_names[j] & Species == i]) > 0){
        amt=as.numeric(spec_abunds[Sample==samp_names[j] & Species == i, value]) #amount in that sample
        if(amt>0){
          kolist=rbind(kolist, data.table(KO=as.character(ref_kos[Species==i,KO]),Abund=ref_kos[Species==i,as.numeric(CopyNum)]*amt/scale_factor, Species=i, Sample=samp_names[j]))
          ##scale by 1000 to make numbers more manageable
        }
      } #else kos_by_sample[[j]]=NULL
    }
  }
  all_kos = kolist[,sum(Abund),by=list(KO,Sample)]
  all_kos = dcast.data.table(all_kos,KO~Sample, value.var="V1")
  return(all_kos)
}


#' Test for enrichment of list of known microbial metabolites in results
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


##Functions for selecting approximate compound identifications based on MetaboSearch output
select_best_id2 = function(met_table2, met_data, net_compounds, final_method = "first"){ ###no retention time, swedish data format
  met_table2[,Mass:=as.character(Mass)]
  met_table2 = melt(met_table2, id.vars = c("Mass", "Delta"))
  subjects = names(met_data)[names(met_data)!="Mass"]
  met_data$Mass = as.character(met_data$Mass)
  met_table2 = met_table2[value %in% net_compounds]
  met_data = met_data[Mass %in% met_table2$Mass]
  met_table2 = met_table2[Mass %in% met_data$Mass,]
  met_table2 = unique(met_table2[order(Mass),])
  unique_ids = sort(unique(met_table2[,value]))
  good_mets = c()
  new_mets = data.table()
  for(j in 1:length(unique_ids)){
    foo = met_table2[value==unique_ids[j]]
    foo_good = foo[which.min(foo$Delta),]
    foo_data = met_data[Mass %in% foo$Mass,]
    #foo_good_close = foo[abs(foo$RetTime - foo_good$RetTime)<0.1,]
    if(dim(foo_data)[1] == 1){
      new_mets = rbind(new_mets, foo_data[,subjects,with=F])
      #figure out way to merge ones that are close in mass
      good_mets = c(good_mets, unique_ids[j])
      met_data = met_data[!Mass %in% foo[,Mass]] #remove these peaks from consideration for future identifications
      met_table2 = met_table2[!Mass %in% foo[,Mass]]
    }
#     if(nrow(foo_data) > 1){
#       test = apply(foo_data[,1:ncol(foo_data)],2,function(x){length(x[x!=0])})
#       if(length(test[test>1])>2){
#         foo_good=foo[which.min(foo$Delta),]
#         foo_data2=foo_data[foo_data$Mass %in% foo_good$Mass,subjects]
#       } else foo_data2=apply(foo_data[,1:12],2,sum)
#       new_mets = rbind(new_mets, foo_data2)
#       good_mets = c(good_mets,unique_ids[j])
#       met_data = met_data[-which(met_data$Mass %in% foo$Mass),] #remove these peaks from consideration for future identifications
#       met_table2 = met_table2[-which(met_table2$Mass %in% foo$Mass),]
#     }
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
  #most_rxns not yet implemented s
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

#Capitalize 1st letter of words
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

#' Defaults for NULL values
#' @export
`%||%` <- function(a, b) if (is.null(a)) b else a


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



