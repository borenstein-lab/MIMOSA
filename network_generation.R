##Generating full network using mapformula, filtering changed values, processing stoichiometry
options(stringsAsFactors = F)
library(data.table)
args = commandArgs(trailingOnly = T)

#reaction_mapformula.lst file from the KEGG ftp server is required to get only reactions annotated in pathways.
#We are just using mapformula plus general reaction-KO mapping - this is basically following the same process as was used to generate the above network
  #mapformula for minpath/pathways
mapformula = fread("reaction_mapformula.lst", colClasses = "character") #get mapformula pathway annotations of reactions
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
write.table(mapformula, file = "reaction_mapformula_split.txt", sep = "\t", quote=F, row.names=F)

###Combine processed mapformula with information on reaction IDs and stoichiometry
#### Option 1: Access most recent KEGG reaction annotations with the KEGGREST API
if(args[1]=="KEGGREST"){
  library(KEGGREST)
  kos_rxns = lapply(unique(mapformula[,Rxn]), keggLink, target="ko")
  kos_to_rxns = data.table(Rxn = gsub("rn:","", unlist(sapply(kos_rxns, function(x){ return(names(x))}))), KO = gsub("ko:","",unlist(kos_rxns)))
  all_kegg = list(KOs = gsub("ko:","",names(keggList("ko")), fixed=T),
                  Reactions = gsub("rn:", "", names(keggList("reaction")), fixed = T))
  all_kegg$Reaction_info = lapply(all_kegg$Reaction, keggGet)
} else {
  ##### Option 2: Read reaction info and KO links from KEGG database file
  kos_to_rxns = fread("KEGG/genes/ko/ko_reaction.list", header=F)
  setnames(kos_to_reactions, c("KO", "Rxn"))
  kos_to_rxns[,KO:=gsub("ko:","",KO)]
  kos_to_rxns[,Rxn:=gsub("rxn:","", Rxn)]
  all_kegg = list(KOs = kos_to_rxns[order(Rxn)][,unique(KO)], Reactions = kos_to_rxns[order(Rxn)][,unique(Rxn)])
  library(readr)
  reaction_info = strsplit(read_file("KEGG/ligand/reaction/reaction"), split = "///\n")[[1]]
  assoc_ids = gsub(".*ENTRY(.*)Reaction.*","\\1",reaction_info)
  all_kegg$Reaction_info = reaction_info[which(assoc_ids %in% all_kegg$Reactions)]
}


mapformula = mapformula[Path!=" 01100"]
mapformula = merge(mapformula, kos_to_rxns, by="Rxn", all.x=T, all.y=F, allow.cartesian=T)
mapformula = mapformula[!is.na(KO)]
setkey(mapformula, NULL)
    rxn_table = unique(mapformula)
    #rxn_table = merge(rxn_table, mapformula, all.x=T, all.y=F, by=c("Rxn","Reac","Prod"), allow.cartesian = T)
    #setkey(rxn_table,NULL)
    #rxn_table = unique(rxn_table)
    #pull out rxn ids from mapformula that include the same reac-prod combo - probably more efficient ways to do this

    #load("KeggReactions.rda")

    #save only reactions that actually involve the correct ko - not sure how these got in here
    #rxn_table[,rxn_id:=rxn_id[sapply(rxn_id, function(y){ KO %in% names(all_kegg$Reaction_info[[match(y, all_kegg$Reactions)]]$ORTHOLOGY)})]]
    #Only include KO-Rxn-Compound connections that are consistent across all files
    rxn_table_sub = unique(rxn_table[,list(Rxn,KO,Reac,Prod)])
    rxn_id_check = sapply(1:length(rxn_table_sub[,Rxn]), function(x){
      ind = which(all_kegg$Reactions==rxn_table_sub[x,Rxn])
      if(length(ind) < 1){ return("noMatch")
      } else if(!(rxn_table_sub[x,KO] %in% names(all_kegg$Reaction_info[[ind]]$ORTHOLOGY))){
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
      eqn = strsplit(rxn_info$EQUATION, split = " <=> ", fixed=T) #get chemical formula for reaction
      comps = unlist(lapply(eqn, strsplit, split = " \\+ "))
      coefs = unique(comps[grepl(" ", comps)])
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
    write.table(rxn_table, file = "ko_rxn_map_all_info_filtered.txt", sep = "\t", quote=F, row.names=F)
