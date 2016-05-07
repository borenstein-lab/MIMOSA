##Generating full network using mapformula, filtering changed values, processing stoichiometry

#     fc = file("network_map_ignoreP01100_rc_to_edges.txt.unix")
#     map = readLines(fc)
#     close(fc)
#     map = strsplit(map, "\t")
#     map = lapply(map, function(x){ return(c(x, rep('',9-length(x))))})
#     map = do.call("rbind",map)
#     map = data.table(map)
#     setnames(map, c("Rxn","R1","R2","R3","R4","R5","R6","R7","R8"))
#     map = melt(map, id.var="Rxn")
#     map = map[value!='']
#     map[,variable:=NULL]
#     maprxns = strsplit(map[,value], split= " pp ")
#     map[,Reac:=sapply(maprxns, function(x){ return(x[1])})]
#     map[,Prod:=sapply(maprxns, function(x){ return(x[2])})]
#     map[,value:=NULL]
#     compounds=unique(c(map[,Reac], map[,Prod]))
#     kos_to_rxns = fread("ko_rxn_ids.txt")
#     map = merge(map, kos_to_rxns, by="Rxn", all.x=T, all.y=F, allow.cartesian = T)
#     map = map[!is.na(KO)]
#     map = map[Reac != Prod]
#     write.table(map, file = "network_map_ignoreP01100_rc_ko_edges_all.txt", sep = "\t", quote=F, row.names=F)
#     
    

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
    #I think it might be best just to use mapformula - I think they are basically the same
    #rxn_table = fread("network_map_ignoreP01100_rc_ko_edges_all.txt")    
    mapformula = fread("reaction_mapformula_split.txt", colClasses = "character")
    mapformula = mapformula[Path!=" 01100"]    
    mapformula = merge(mapformula, kos_to_rxns, by="Rxn", all.x=T, all.y=F, allow.cartesian=T)
    mapformula = mapformula[!is.na(KO)]
    setkey(mapformula, NULL)
    rxn_table = unique(mapformula)
    #rxn_table = merge(rxn_table, mapformula, all.x=T, all.y=F, by=c("Rxn","Reac","Prod"), allow.cartesian = T) 
    #setkey(rxn_table,NULL)
    #rxn_table = unique(rxn_table)
    #pull out rxn ids from mapformula that include the same reac-prod combo - probably more efficient ways to do this    
    
    load("KeggReactions.rda")
    
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
