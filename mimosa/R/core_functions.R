#core_functions.R
#workflow of functions for MIMOSA analysis


get_kegg_reaction_links = function(rxns){
  foo = vector("list", length(rxns))
  for(j in 1:length(rxns)){ ##Using lapply here causes weird timeouts with curl
    foo[[j]] = tryCatch(KEGGREST::keggLink(rxns[j], target="ko"))
  }
  return(foo)
}

#' Get basic KEGG reaction and KO info in a unified format
#'
#' @import data.table
#' @param kos_to_rxns_method EITHER "KEGGREST" indicating to use the KEGGREST API to link KOs and reactions and get reaction info, OR a file path to the KEGG file genes/ko/ko_reaction.list
#' @param reaction_info_file If kos_to_rxns_method is a file path, additional file containing full reaction info from the KEGG database
#' @param save_out whether to save output as an Rdata file named "KeggReactions.rda"
#' @param kolist Optionally, a vector of KO IDs. Will create an all_kegg object containing information only on those KOs and the reactions linked to them.
#' @return An R object with 4 components: a list of KOs, a list of Reaction IDs, a list of associated reaction info, and a table linking KOs to reactions
#' @examples
#' get_kegg_reaction_info("KEGGREST")
#' get_kegg_reaction_info("KEGG/genes/ko/ko_reaction.list", "KEGG/ligand/reaction/reaction")
#' @export
get_kegg_reaction_info = function(kos_to_rxns_method, reaction_info_file = "", save_out = T, kolist = ""){
#### Option 1: Access most recent KEGG reaction annotations with the KEGGREST API
  if(kos_to_rxns_method=="KEGGREST"){
      all_kegg = list(KOs = gsub("ko:","",names(KEGGREST::keggList("ko")), fixed=T),
                      Reactions = gsub("rn:", "", names(KEGGREST::keggList("reaction")), fixed = T))
    kos_to_rxns = get_kegg_reaction_links(all_kegg$Reactions)
    all_kegg$Reactions = all_kegg$Reactions[sapply(kos_to_rxns, length) != 0]
    kos_to_rxns = kos_to_rxns[which(sapply(kos_to_rxns, length) != 0)]
    kos_to_rxns = data.table(Rxn = gsub("rn:","", unlist(sapply(kos_to_rxns, function(x){ return(names(x))})), fixed = T), KO = gsub("ko:","",unlist(kos_to_rxns), fixed = T))
    if(!identical(kolist, "")){ #if we only want specific KOs)
      kos_to_rxns = kos_to_rxns[KO %in% kolist]
      all_kegg$KOs = all_kegg$KOs[all_kegg$KOs %in% kolist]
      all_kegg$Reactions = all_kegg$Reactions[all_kegg$Reactions %in% kos_to_rxns[,Rxn]]
    }
    all_kegg$Reaction_info = lapply(all_kegg$Reactions, function(x){
      return(KEGGREST::keggGet(x)[[1]]) }) #This is slow!
    cat("Done downloading reaction info!\n")
  } else {
  ##### Option 2: Read reaction info and KO links from KEGG database file
    kos_to_rxns = fread(kos_to_rxns_method, header=F)
    setnames(kos_to_rxns, c("KO", "Rxn"))
    kos_to_rxns[,KO:=gsub("ko:","",KO)]
    kos_to_rxns[,Rxn:=gsub("rn:","", Rxn)]
    all_kegg = list(KOs = kos_to_rxns[,sort(unique(KO))], Reactions = kos_to_rxns[order(Rxn)][,unique(Rxn)])
    if(!requireNamespace("readr", quietly = T)){
      stop("readr required to process KEGG reaction file, please install", call. = F)
    }
    reaction_info = strsplit(readr::read_file(reaction_info_file), split = "///\n", fixed = T)[[1]]
    assoc_ids = trimws(gsub(".*ENTRY(.*)Reaction.*","\\1",reaction_info))
    all_kegg$Reaction_info = reaction_info[match(all_kegg$Reactions, assoc_ids)]
    all_kegg$Reaction_info = lapply(all_kegg$Reaction_info, function(x){
      x = strsplit( gsub("(\n[A-Z])","~\\1",x), "~\n", fixed = T)[[1]]
      new_ob = list()
      for(j in 1:length(x)){
        foo = strsplit(x[j], split = "  +")[[1]]
        new_ob[[foo[1]]] = foo[2:length(foo)]
      }
      return(new_ob)
    })
  }
  all_kegg$kos_to_rxns = kos_to_rxns
  if(save_out) save(all_kegg, file = "KeggReactions.rda")
  return(all_kegg)
}

#' Generate community metabolic network template from the KEGG database, using reaction_mapformula.lst and
#'
#' @import data.table
#' @param mapformula_file file path to "reaction_mapformula.lst" from the KEGG database, specifying reaction IDs annotated in pathways
#' @param all_kegg Output of get_kegg_reaction_info
#' @return A table specifying a community metabolic network based on KEGG pathways with the following columns: Rxn (reaction ID),	KO (Gene ID),	Reac (Reactant ID),	Prod (Product ID),	Path (Pathway ID),	ReacProd (Original reaction specification),	stoichReac (Reactant stoichiometry coefficient),	stoichProd (Product stoichiometry coefficient)
#' @examples
#' generate_network_template_kegg("KEGG/ligand/reaction/reaction_mapformula.lst", all_kegg)
#' @export
generate_network_template_kegg = function(mapformula_file, all_kegg, save_out = F){
  mapformula = fread(mapformula_file, colClasses = "character") #get mapformula pathway annotations of reactions
  setnames(mapformula, c("Rxn","Path","ReacProd"))
  #Process Reacs and Prods, flip reversible reactions, etc
  mapformula[,Reac:=lapply(ReacProd,function(x){ return(unlist(strsplit(x,"="))[1])})]
  mapformula[,Reac:=lapply(Reac, strsplit, " ")]
  mapformula[,Prod:=lapply(ReacProd,function(x){ return(unlist(strsplit(x,"="))[2])})]
  mapformula[,Prod:=lapply(Prod, strsplit, " ")]
  mapformula[,Reac:=sapply(Reac, unlist)]
  mapformula[,Prod:=sapply(Prod, unlist)]
  mapformula[,Reac:=sapply(Reac, function(x){ return(x[grepl("[A-Za-z]",x)])})]
  mapformula[,Prod:=sapply(Prod, function(x){ return(x[grepl("[A-Za-z]",x)])})]
  for(j in 1:length(mapformula[,Rxn])){
    submap = mapformula[j]
    if(grepl("<=>",submap[,ReacProd])){ #add other way if flipped
      setnames(submap, c("Reac","Prod"),c("Prod","Reac"))
      submap = submap[,list(Rxn,Path,ReacProd,Reac,Prod)]
      mapformula = rbind(mapformula, submap)
    }
    if(grepl("<= ",submap[,ReacProd])){ #flip if rxn backward
      setnames(submap, c("Reac","Prod"),c("Prod","Reac"))
      submap = submap[,list(Rxn,Path,ReacProd,Reac,Prod)]
      mapformula[j] = submap
    }
  }
  new_mapformula = data.table(Rxn=rep('',0),Path=rep('',0), ReacProd=rep('',0), Reac=rep('',0),Prod=rep('',0))
  for(j in 1:length(mapformula[,Rxn])){ #add extra rows for multiple reacs and prods
    reacs = mapformula[j,Reac][[1]]
    prods = mapformula[j,Prod][[1]]
    submap = mapformula[j]
    for(k in 1:length(reacs)){
      for(m in 1:length(prods)){
        submap[,Reac:=reacs[k]]
        submap[,Prod:=prods[m]]
        new_mapformula = rbind(new_mapformula, submap)
      }
    }
  }
  mapformula = new_mapformula
  mapformula = mapformula[Path!=" 01100"]
  mapformula = merge(mapformula, all_kegg$kos_to_rxns, by="Rxn", all.x=T, all.y=F, allow.cartesian=T)
  mapformula = mapformula[!is.na(KO)]
  setkey(mapformula, NULL)
  rxn_table = unique(mapformula)
  rxn_table_sub = unique(rxn_table[,list(Rxn,KO,Reac,Prod)])
  rxn_id_check = sapply(1:length(rxn_table_sub[,Rxn]), function(x){
    ind = which(all_kegg$Reactions==rxn_table_sub[x,Rxn])
    if(length(ind) < 1){ return("noMatch")
    } else if(!(rxn_table_sub[x,KO] %in% names(all_kegg$Reaction_info[[ind]]$ORTHOLOGY)|rxn_table_sub[x,KO] %in% all_kegg$Reaction_info[[ind]]$ORTHOLOGY)){
      return("noKOmatch")
    }
    else if(!(grepl(rxn_table_sub[x,Reac], all_kegg$Reaction_info[[ind]]$EQUATION) & grepl(rxn_table_sub[x,Prod], all_kegg$Reaction_info[[ind]]$EQUATION))){
      return("noCmpdmatch")
    }
    else return("good")
  })
  rxn_table_sub = rxn_table_sub[rxn_id_check=="good"]
  rxn_table = merge(rxn_table, rxn_table_sub, by=c("Rxn","KO","Reac","Prod"), all.x=F)

  #get stoichiometry too
  goodrxns = match(rxn_table[,Rxn], all_kegg$Reactions)
  stoichReac = rep(1, length(rxn_table[,KO]))
  stoichProd = rep(1, length(rxn_table[,KO]))
  #parse KEGG Reaction for stoichiometry
  for(j in 1:length(goodrxns)){
    rxn_info = all_kegg$Reaction_info[[goodrxns[j]]]
    if(!is.null(rxn_info$EQUATION)){
      eqn = strsplit(rxn_info$EQUATION, split = " <=> ", fixed=T) #get chemical formula for reaction
      comps = unlist(lapply(eqn, strsplit, split = " \\+ "))
      coefs = unique(comps[grepl(" ", comps)])
    } else {
      coefs = c()
    }
    if(length(coefs) > 0){
      coefs_split = strsplit(coefs, split = " ",fixed=T)
      comp = sapply(coefs_split, function(x){ return(x[2])})
      reac_coef = unlist(coefs_split[comp==rxn_table[j,Reac]])[1]
      if(length(unlist(reac_coef)) > 0){
        if(!is.na(as.numeric(reac_coef))) stoichReac[j] = as.numeric(reac_coef)
      }
      prod_coef = unlist(coefs_split[comp==rxn_table[j,Prod]])[1]
      if(length(unlist(prod_coef)) > 0)
        if(!is.na(as.numeric(prod_coef))) stoichProd[j] = as.numeric(prod_coef)
    }
  }
  rxn_table[,stoichReac:=stoichReac]
  rxn_table[,stoichProd:=stoichProd]
  rxn_table[,Path:=gsub(" ","",Path)]
  if(save_out) write.table(rxn_table, file = "mapformula_all_info.txt", quote=F, row.names = F, sep = "\t")
  return(rxn_table)
}

#' Create a community metabolic network model using a few different methods.
#'
#' @import data.table
#' @param kos Genes to include in network (KEGG Orthology IDs)
#' @param keggSource source of network information, currently can be "KeggTemplate", "loadNet", or "metacyc"
#' @param degree_filter Compounds connected to this number of KOs or more will be filtered from the network
#' @param minpath_file file of minimal pathways to use for core network
#' @param normalize Whether to normalize output matrix to show relative impacts of each gene on each compound
#' @param rxn_table If keggSource is "KeggTemplate", must supply a data.table with the format of
#' @param networkFile If keggSource is "loadNet", file that template network should be loaded from, should have same format as output of generate_network_template_kegg
#' @return List containing 3 different versions of the same network: an adjacency matrix, an adjacency matrix that differentiates between genes with a neutral effect on a compound and no effect (0 vs NA), and an edge list.
#' @examples
#' generate_genomic_network(kos, "KeggTemplate", degree_filter = 30)
#' @export
generate_genomic_network = function(kos, keggSource = "KeggTemplate", degree_filter = 0,
                                    minpath_file = '', normalize = T, rxn_table = "", networkFile=""){
  #keggSource must be either "labKegg", or "loadNet" for if the network has already been generated, or "metacyc"
  #degree filter says whether to filter hub compounds with a degree greater than this value
  if(!keggSource %in% c("loadNet", "KeggTemplate", "metacyc")){
    stop("Invalid network method selected")
    } else if(keggSource == "loadNet"){
    load(networkFile)
    return(allnet)
  } else if(keggSource == "KeggTemplate") { #load network template and just grab subset
    if(identical(rxn_table, "")){
      stop("Must supply a community network reaction table, for example the output of generate_network_template_kegg")
    }
    #rxn_table = fread("ko_rxn_map_all_info_filtered.txt", colClasses = c(rep("character",6), rep("numeric",2)))
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
    goodrxns = sapply(all_kegg$Reaction_info, function(x){ any(names(x$ORTHOLOGY) %in% kos) | any(x$ORTHOLOGY %in% kos)})
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
    rxn_ids = rxn_ids[rxn_ids %in% metacyc_rxns[,RID]]
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

    ##get stoichiometry and directionality from metacyc
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
#' @param minpath_file Optional file of pathways to filter network to
#' @param cor_method Either "spearman" or "pearson", default is Spearman
#' @param net_file file containing network template to use
#' @param nperm Number of permutations for Mantel test, default is 20000
#' @param nonzero_filter Minimum number of samples required to have nonzero concentrations and nonzero metabolic potential scores in order for metabolite to be evaluated, default is 3
#' @return No return, writes output to file
#' @examples
#'
#' @export
run_all_metabolites = function(genes, mets, file_prefix = 'net1', correction = "fdr", cutoff = 0.1,
                               net_method = "load", rxn_table_source = "", id_met = F, met_id_file = '',
                               degree_filter = 0, minpath_file = '', cor_method = "spearman",
                               net_file = "", nperm = 20000, nonzero_filter=3){
  #must be data.tables keyed by KOs/KEGG IDs in first column, all other columns are subject IDs
  #correction must be either "bonferroni" or "fdr", cutoff is q value cutoff
  #id_mets specifies whether to use network for improved metabolite identification (i.e. for Braun datasets)
  #if true, text file with list of putative met ids from Metabosearch must be included
  kos = genes[,KO]
  nsamp = dim(genes)[2] - 1
  subjects = intersect(names(genes), names(mets))
  ko_net = generate_genomic_network(kos, keggSource = net_method, degree_filter = degree_filter, minpath_file = minpath_file, rxn_table = rxn_table_source, networkFile = net_file)
  ko_network_mat = ko_net[[1]]
  stoich_mat = ko_net[[2]]
  ko_net_table = ko_net[[3]]
  if(id_met){
    cat("Selecting best metabolite identifications\n")
    net_compounds = row.names(ko_net[[1]])
    met_id_table = fread(met_id_file, sep = "\t", header=T)
    met_id_table = met_id_table[met_id_table$Delta!="-"] #get rid of mets with no possible KEGG ID
    #met_id_list = strsplit(unlist(met_id_list), split = " ")
    #if(grepl("swedish",file_prefix)){
    new_mets = select_best_id2(met_id_table, mets, net_compounds)
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
    if(length(met_vector[met_vector!=0]) < nonzero_filter | length(prmt_vector[prmt_vector!=0]) < nonzero_filter){
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
#' @param minpath_file Optional file of pathways to filter network to
#' @param net_file file containing network template to use
#' @param quant Quantile above which a metabolite is "elevated", default is 0.5
#' @param plot_rank Whether to generate plots of ranks of metabolite concentrations and scores, default is F
#' @param plot_continuous Whether to generate plots of metabolite concentrations and scores, default is F
#' @param nonzero_filter Minimum number of samples required to have nonzero concentrations and nonzero metabolic potential scores in order for metabolite to be evaluated, default is 3
#' @return No return, writes output to file
#' @examples
#'
#' @export
run_all_metabolites2 = function(genes, mets, file_prefix = 'net1', correction = "fdr", cutoff = 0.1,
                               net_method = "load", rxn_table_source = "", id_met = F, met_id_file = '',
                               degree_filter = 0, minpath_file = '',
                               net_file = "", quant = 0.5, plot_rank = F, plot_continuous=F, nonzero_filter=3, rel_abund_mets = F){
  #must be data.tables keyed by KOs/KEGG IDs in first column, all other columns are subject IDs
  #correction must be either "bonferroni" or "fdr", cutoff is q value cutoff
  #id_mets specifies whether to use network for improved metabolite identification (i.e. for Braun datasets)
  #if true, text file with list of putative met ids from Metabosearch must be included
  #can predict above/below median or other threshold with quant argument
  kos = genes[,KO]
  nsamp = dim(genes)[2] - 1
  subjects = intersect(names(genes), names(mets))
  ko_net = generate_genomic_network(kos, keggSource = net_method, degree_filter = degree_filter, minpath_file = minpath_file,
                                    rxn_table = rxn_table_source, networkFile = net_file)
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
    if(length(met1[met1!=0]) < nonzero_filter | length(prmts[prmts!=0]) < nonzero_filter){
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
#' @param j index (from lapply usually)
#' @param prmts_sub_good CMP scores for metabolites to be analyzed
#' @param all_rxns list of reaction tables for each compound
#' @param subjects vector of subject names
#' @param norm_kos gene abundance matrix
#' @param ko_net full network
#' @return list of contributors and correlations
#' @examples
#' lapply(1:length(good_mets), prmt_contributions, prmts_sub_good = prmts_sub_good, all_rxns = all_rxns[[j]],
#' subjects=subjects, norm_kos = norm_kos, ko_net = ko_net)
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


#' Randomize a metabolic network edge list by randomly sampling 2 edges and if it works, switching products
#'
#' @import data.table
#' @param netw Network edge list, format of the 3rd output item of generate_genomic_network
#' @param n_reps Number of successful edge switches to perform (default = 5000)
#' @return Network edge list in the same format with the same compounds and degree distribution, but with randomized edges
#' @examples
#' randomize_net(edge_list, 3000)
#' @export
randomize_net=function(netw, n_reps = 5000){
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

#' Run one iteration of MIMOSA analysis with a randomized metabolic network. Run this function many times to generate a null distribution of metabolite comparisons.
#'
#' @import data.table
#' @param out_file File path to .Rda file of MIMOSA output from run_all_metabolites
#' @param id_num ID number of shuffle (if running many iterations)
#' @param n_iter Number of edge switches to generate randomized network, default is 5000
#' @param nonzero_filter Minimum number of samples required to have nonzero concentrations and nonzero metabolic potential scores in order for metabolite to be evaluated, default is 3
#' @param qval_thresholds 1 or more significance thresholds above which to count the number of metabolites
#' @return None, writes two tables to file - one of the total counts of metabolites meeting each threshold, and one of the results for each metabolite
#' @examples
#' run_shuffle("Dataset2_bv_out.rda", id_num = 1)
#' @export
#'
run_shuffle = function(out_file, id_num = 1, n_iter = 5000, nonzero_filter = 3, qval_thresholds = c(0.1, 0.01)){
  load(out_file)
  rxn_table = ko_net[[3]]
  random_net = randomize_net(rxn_table, n_iter)
  net_mats = make_network_matrix(random_net)
  prmts = get_prmt_scores(net_mats[[1]], norm_kos)

  metIDs = mets[,KEGG]
  shared_mets = metIDs[metIDs %in% row.names(net_mats[[1]])]
  mets_shared_only = mets[shared_mets]
  prmts_shared_only = prmts[shared_mets]

  good_data = which(apply(prmts_shared_only, 1, function(x){ length(x[as.numeric(x)!=0]) >= nonzero_filter }) & apply(mets_shared_only, 1, function(x){ length(x[as.numeric(x)!=0]) >= nonzero_filter }))
  mets_shared_only = mets_shared_only[good_data]
  prmts_shared_only = prmts_shared_only[good_data]
  shared_mets = shared_mets[good_data]

  compare_pos = sapply(1:length(shared_mets), function(x){
    compare_met(shared_mets[x], shared_mets[x], mets_shared_only, prmts_shared_only, posneg = "pos", nperm = 10000)
  })
  compare_neg = sapply(1:length(shared_mets), function(x){
    compare_met(shared_mets[x], shared_mets[x], mets_shared_only, prmts_shared_only, posneg = "neg", nperm = 10000)
  })

  corrected_pos = correct(compare_pos, method="fdr")
  corrected_neg = correct(compare_neg, method = "fdr")
  total_pos = c()
  total_neg = c()
  for(j in 1:length(qval_thresholds)){
    total_pos[j] = length(corrected_pos[corrected_pos < qval_thresholds[j] & compare_pos < qval_thresholds[j] & !is.na(corrected_pos)])
    total_neg[j] = length(corrected_neg[corrected_neg < qval_thresholds[j] & compare_neg < qval_thresholds[j] & !is.na(corrected_neg)])
  }

  met_table = data.table(compound=shared_mets, Pos=corrected_pos, Neg=corrected_neg, Iter=rep(id_num, length(shared_mets)))

  file_prefix = paste0(gsub("_out.rda", "", out_file), id_num)

  write.table(data.table(t(total_pos), t(total_neg)), file = paste0(file_prefix, "_shuff_network.txt"), sep = "\t", quote=F, row.names=F, col.names=F, header = F)
  write.table(met_table, file = paste0(file_prefix,"_shuff_network_mets.txt"), sep = "\t", quote=F, row.names=F, col.names=F)
}

