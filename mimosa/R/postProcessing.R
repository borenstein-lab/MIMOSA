###Post-processing functions
###For summarizing and visualizing output of the MIMOSA model


#' Read output of runMimosa.R from files to make a table of all relevant MIMOSA results for metabolites identified as consistent
#'
#' @import data.table
#' @param prefix Naming prefix used by runMimosa.R
#' @param contribs_file File of species-specific gene abundances used for calculating species contributions
#' @return A table where every row represents genes and species that were identified as potential important contributors for a consistent metabolite
#' get_all_consistent_metabolite_info("run1", "filePath/otu_gene_contributions.txt")
#' @export
get_all_consistent_metabolite_info = function(prefix, contribs_file){
  consistent_node_data = fread(paste0(prefix,"_signifPos.txt"))
  setnames(consistent_node_data, "Correlation", "CMP_Correlation")
  #consistent_node_data[,CorrN:=NULL]
  consistent_node_data[,Metabolite:=met_names(compound)] #This will tell you metabolite names for each ID
  species_data = fread(paste0(prefix,"_specContrib.txt"))
  consistent_node_data = merge(consistent_node_data, species_data[Pass==1 & !is.na(Pass)], by="compound", all.x = T, all.y = F) # Integrate species data into same table
  setnames(consistent_node_data, "Cor", "SpeciesCor")
  ko_summary = fread(paste0(prefix, "_geneContribCompoundSummary.txt"))
  consistent_node_data = merge(consistent_node_data, ko_summary, by="compound", all.x = T, all.y = F)

  #Process gene contributions analysis to only include the "important" genes:
  geneContribs = fread(paste0(prefix,"_geneContribAnalysis.txt"))
  geneContribs = geneContribs[Cor < 0.5|is.na(Cor)] #This will retain only potential "key contributor" genes using the same threshold heuristic we used in our paper

  ##You can also merge this into the main table above, but have to first make it species-specific:
  species_ref = unique(fread(contribs_file)[,list(Gene,OTU)])
  setnames(species_ref, c("KO", "Species"))
  geneContribs = merge(geneContribs, species_ref, by = "KO", all.x = T, all.y = F)

  consistent_node_data = merge(consistent_node_data, geneContribs, by = c("compound", "Species"), all.x = T, all.y = F, allow.cartesian = T)
  setnames(consistent_node_data, "Cor", "GeneCor")

  return(consistent_node_data)
}

#' Add summary columns to metabolite 'node' file of MIMOSA results
#'
#' @import data.table
#' @param node_data Table from "_nodes.txt" file produced by runMimosa.R
#' @param thresholds Significance thresholds (q-value and p-value) to be used for classifying metabolites as consistent or contrasting. By default this function classifies based on 2 different thresholds
#' @param ko_summary Table of gene contributor information produced by MIMOSA (geneContribCompoundSummary.txt))
#' @param path_key Table of metabolite categories
#' @return Processed table of merged results for each metabolite
#' process_metabolite_summary(node_data, list(c(0.01, 0.01), c(0.1, 0.1)))
#' @export
process_metabolite_summary = function(node_data, thresholds, ko_summary){
  threshold1 = thresholds[[1]]
  threshold2 = thresholds[[2]]
  node_data[,Consistent:=ifelse(QValPos < threshold1[1] & PValPos < threshold1[2], 1,0)]
  node_data[,Contrasting:=ifelse(QValNeg < threshold1[1] & PValNeg < threshold1[2], 1,0)]
  node_data[,PredictionType:=ifelse(Consistent==1, "Consistent", "Inconsistent")]
  node_data[,PredictionType:=ifelse(Contrasting==1, "Contrasting", PredictionType)]
  node_data[,PredictionType2:=ifelse(QValPos < threshold2[1] & PValPos < threshold2[2], "Consistent", "Inconsistent")]
  node_data[,PredictionType2:=ifelse(QValNeg < threshold2[1] & PValNeg < threshold2[2], "Contrasting", PredictionType2)]
  node_data[,PredictionType2:=factor(PredictionType2, levels = c("Inconsistent", "Contrasting", "Consistent"))]
  node_data[,PredictionType:=factor(PredictionType, levels = c("Inconsistent", "Contrasting", "Consistent"))]
  #Add gene summary data
  node_data = merge(node_data, ko_summary, by="compound", all.x = T, all.y = F)
  node_data[,PrimaryMake2:=ifelse(PrimaryMake==1, "Synthesis", "Combination")]
  node_data[,PrimaryMake2:=ifelse(PrimaryMake==-1, "Degradation", PrimaryMake2)]
  node_data[,PrimaryMake2:=factor(PrimaryMake2, levels = c("Synthesis", "Combination", "Degradation"))]
  #Add metabolite category information
  node_data = merge(node_data, path_key, by ="compound", all.x = T, all.y = F)
  node_data[CompoundName=="D-Galactose", SuperPath:="Carbohydrate"]
  node_data[CompoundName=="Serotonin", SuperPath:="Amino acid"]
  node_data[CompoundName=="alpha,alpha-Trehalose", SuperPath:="Carbohydrate"]
  met_order = node_data[order(Correlation, decreasing = F),Metabolite]
  node_data[,Metabolite:=factor(Metabolite, levels = met_order)]
  return(node_data)
}

#' Do basic enrichment testing of between-group metabolite differences and add results to the MIMOSA summary table
#'
#' @import data.table
#' @param node_data Table from "_nodes.txt" file produced by runMimosa.R
#' @param mets Table of metabolite abundance data
#' @param metadata Table of samples and their grouping by a binary variable
#' @param metadata_var Column name in metadata table for the binary grouping variable to be used for testing
#' @return Processed table including columns describing binary enrichment status and difference in mean
#' add_metadata_to_metabolite_summary(node_data, mets, metadata, "DiseaseStatus")
#' @export
add_metadata_to_metabolite_summary = function(node_data, mets, metadata, metadata_var){
  mets = merge(mets, metadata[,c("Sample", metadata_var), with=F], all.x = T)
  bad_mets = mets[!is.na(value), length(Sample),by=c("KEGG",metadata_var)]
  missing_mets = bad_mets[,length(unique(get(metadata_var))), by=KEGG]
  bad_mets = c(bad_mets[V1 < 3, unique(KEGG)], missing_mets[V1 < 2, KEGG])
  mets_summary = mets[!is.na(value) & !KEGG %in% bad_mets,list(mean(value[get(metadata_var)==1])-mean(value[get(metadata_var)==0]),wilcox.test(value[get(metadata_var)==1],value[get(metadata_var)==0])$p.value),by=KEGG]
  setnames(mets_summary, c("V1", "V2"), c("MetDiff", "MetPVal"))
  node_data = merge(node_data, mets_summary, by.x = "compound", by.y = "KEGG", all.x = T, all.y = F)
  node_data[,Enriched:=ifelse(MetDiff > 0 & MetPVal < 0.1,1,0)]
  node_data[,Depleted:=ifelse(MetDiff < 0 & MetPVal < 0.1,1,0)]
  node_data[,Status:=ifelse(Enriched, "Enriched", "None")]
  node_data[,Status:=ifelse(Depleted, "Depleted", Status)]
  node_data[,Status:=factor(Status, levels = c("Enriched", "None", "Depleted"))]
  return(node_data)
}

#' Make table summarizing species contributions for each analyzed metabolite
#'
#' @import data.table
#' @param species_data Table from MIMOSA species contributor analysis ("_specContrib.txt")
#' @param node_data Table of core metabolite results
#' @param tax_ref Table of taxonomy assignments for each species
#' @return Table describing species contributors for each metabolite along with its overall MIMOSA results
#' make_metabolite_species_contribution_table(species_data, node_data, tax_ref, path_key)
#' @export
make_metabolite_species_contribution_table = function(species_data, node_data, tax_ref, path_key){
  node_data_spec = merge(node_data, species_data[Pass==1 & !is.na(Pass)], by="compound", all.x = T, all.y = F) # Integrate species data into same table
  setnames(node_data_spec, "Cor", "SpeciesCor")
  node_data_spec = merge(node_data_spec, tax_ref, by = "Species", all.x = T, all.y = F)
  node_data_spec[,GenusShort:=gsub("k__.*f__","", Genus)]
  node_data_spec[GenusShort=="; g__", GenusShort:=gsub("k__.*p__", "", Genus)]
  node_data_spec = node_data_spec[!is.na(Species)]
  met_order = node_data[order(Correlation, decreasing = F),Metabolite]
  node_data_spec[,Metabolite:=factor(Metabolite, levels = met_order)]
  return(node_data_spec)
}


#' Process table of gene and reaction contributors
#'
#' @import data.table
#' @param geneContribs Table from MIMOSA gene/rxn contributor analysis ("_geneContribAnalysis.txt")
#' @param node_data Table of core metabolite results
#' @return Merged table of MIMOSA results and gene and reaction contributors
#' process_gene_contribs(gene_contribs, node_data)
#' @export
process_gene_contribs = function(geneContribs, node_data){
  geneContribs = merge(geneContribs, node_data, by = "compound", all.x = T)
  geneContribs[,Contrib:=ifelse(is.na(Cor)|Cor < 0.5, 1, 0)]
  noContribs = geneContribs[,sum(Contrib),by=compound][V1==0,compound]
  geneContribs[compound %in% noContribs, Contrib:=2] #Mark contributors of compounds with no notable ones
  return(geneContribs)
}

#' Make table of species-specific gene and reaction contributors identified by MIMOSA
#'
#' @import data.table
#' @param geneContribs Table from MIMOSA gene/rxn contributor analysis ("_geneContribAnalysis.txt")
#' @param species_ref Table of genome content as produced by PICRUSt
#' @param node_data_spec Table of core metabolite results and their species contributors
#' @return Table of species-specific reaction contributors and overall MIMOSA results for each analyzed metabolite
#' make_gene_species_contributor_table(gene_contribs, species_ref, node_data_species)
#' @export
make_gene_species_contributor_table = function(geneContribs, species_ref, node_data_spec){
  geneContribs2 = merge(geneContribs, species_ref, by = "KO", all.x = T, all.y = F) #Get all species that have that gene
  gene_species = merge(node_data_spec, geneContribs2, by = intersect(names(node_data_spec), names(geneContribs2)), all.x = T, all.y = F, allow.cartesian = T)
  setnames(gene_species, "Cor", "GeneCor")
  return(gene_species)
}

#' Make heatmap plot of metabolite concentrations across samples for analyzed compounds
#'
#' @import data.table
#' @import ggplot2
#' @param mets Table of metabolite abundances
#' @param met_list List of metabolite IDs to include in the plot
#' @return plot of metabolite concentrations across samples
#' plot_metabolite_concentrations(mets)
#' @export
plot_metabolite_concentrations = function(mets, met_list){
  mets[,Sample:=factor(Sample, levels = sort(unique(as.character(Sample))))]
  return(ggplot(mets[KEGG %in% met_list], aes(x=Sample, y = ifelse(is.na(met_names(KEGG)), KEGG, met_names(KEGG)), fill = value)) + geom_tile() + theme(axis.ticks = element_blank(), legend.position = "bottom", axis.text.x = element_text(angle =90, hjust=0), axis.text.y = element_text(size=6)) + ylab("")+ scale_fill_gradient(low = brewer.pal(9, "Blues")[1], high = brewer.pal(9, "Blues")[6]))
}

#' Make bar plot of MIMOSA outcome counts in different metabolite categories
#'
#' @import data.table
#' @import ggplot2
#' @import cowplot
#' @param node_data Processed core MIMOSA results table
#' @param variable Grouping variable for metabolites (must be a column name in node_data)
#' @param threshold Whether to use the lower or higher significance threshold
#' @param prediction_colors List of 3 colors to use as color scheme for Consistent, Contrasting, & Inconsistent metabolites
#' @return plot object of MIMOSA outcome counts across specified metabolite categories
#' plot_metabolite_counts(node_data, "Category")
#' @export
plot_metabolite_counts = function(node_data, variable, threshold = "low", prediction_colors = c("#1B9E77","lightgrey","#F46D43")){
  fill_var = ifelse(threshold=="low", "PredictionType2", "PredictionType")
  plot1 = ggplot(node_data, aes_string(x=variable, fill = fill_var)) + geom_bar(stat="count") + scale_fill_manual(values = c("white", prediction_colors[c(3,1)])) + theme_cowplot() + ylab("Number of metabolites") + theme(legend.title = element_blank(), axis.text.x = element_text(angle=90, hjust=1), axis.ticks = element_blank()) + xlab("") #+ facet_wrap(~Dataset)
  return(plot1)
}

#' Return summarized and concise data on taxonomic contributors for consistent metabolites
#'
#' @import data.table
#' @param node_data_spec Table of species contributors and core MIMOSA results
#' @return Table of summarized and rounded data on taxonomic contributors for consistent metabolites
#' @export
print_taxa_contributor_table = function(node_data_spec){
  if("Status" %in% names(node_data_spec)){ #Allow for not specifying metadata variable
    return(unique(node_data_spec[PredictionType2=="Consistent", list(Correlation = round(Correlation, 3), QValPos = round(QValPos,3), PredType = PredictionType2, Status, Primary = PrimaryMake2, TaxaContrib = paste0(unique(GenusShort), collapse = " ")), by=Metabolite])[order(Correlation, decreasing = T)])
  } else {
    return(unique(node_data_spec[PredictionType2=="Consistent", list(Correlation = round(Correlation, 3), QValPos = round(QValPos,3), PredType = PredictionType2, Primary = PrimaryMake2, TaxaContrib = paste0(unique(GenusShort), collapse = " ")), by=Metabolite])[order(Correlation, decreasing = T)])
  }
}

#' Make plot of individual metabolites and their taxonomic contributors
#'
#' @import data.table
#' @import ggplot2
#' @import cowplot
#' @param node_data Table of core MIMOSA results
#' @param node_data_spec Table of MIMOSA taxonomic contributors
#' @param prediction_colors List of 3 colors for generating overall correlation color scale
#' @return Plot object showing overall MIMOSA correlation and individual taxonomic contributors for each metabolite
#' @export
taxonomic_contributor_heatmap_plot_grid = function(node_data, node_data_spec, prediction_colors = c("#1B9E77","lightgrey","#F46D43")){
  pred_grid1 = ggplot(node_data, aes(x=1, y = Metabolite, fill = Correlation)) + geom_tile(col="black") + theme(axis.line=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), legend.position="bottom", panel.border=element_blank(), axis.text=element_blank()) +scale_fill_gradientn(colours=rev(prediction_colors), limits = c(-1*max(abs(node_data[,Correlation]), na.rm = T),max(abs(node_data[,Correlation]), na.rm = T)))+ scale_x_discrete( expand = c(0, 0)) + guides(size=F, col=F, fill = guide_colourbar(title="Prediction Level", title.position="top")) #+ facet_wrap(~Dataset, nrow = 2)
  spec_order = node_data_spec[,length(Metabolite), by=GenusShort][order(V1, decreasing=T)][,GenusShort]

  #Count number of passing species
  spec_data_counts = node_data_spec[,length(unique(Species)),by=list(Metabolite,GenusShort,PredictionType,PredictionType2)]

  tot_species = unique(node_data_spec[,list(Species,Metabolite)])[,length(unique(Species)), by=Metabolite] #what's the share from this taxon for this metabolite
  spec_data_counts = merge(spec_data_counts, tot_species, by="Metabolite")
  spec_data_counts[,OTUShare:=V1.x/V1.y]
  spec_data_counts = melt(dcast(spec_data_counts, Metabolite+PredictionType+PredictionType2~GenusShort,value.var="OTUShare", fill = 0), id.vars = c("Metabolite", "PredictionType", "PredictionType2"), variable.name = "GenusShort") #fill in missing 0s
  spec_data_counts[,GenusShort:=factor(GenusShort, levels = spec_order)]
  spec_grid = ggplot(spec_data_counts, aes(y=Metabolite, x = GenusShort, fill = value)) + geom_tile(col="grey") + scale_fill_gradient(low = "white", high = brewer.pal(9, "Blues")[9]) + theme(axis.ticks= element_blank(), axis.title = element_blank(), axis.text.x = element_text(angle = 90, hjust=1, size=7), axis.text.y = element_text(size=7), plot.background = element_rect(color = "grey"), panel.grid = element_line(color="gray"), panel.ontop = T) + scale_y_discrete(drop = T) + scale_x_discrete(drop = T) + guides(fill = guide_legend(title = "Share of Contributing OTUs")) #  facet_wrap(~Dataset+PredictionType2, scales = "free_x")

  plot_obj = plot_grid(pred_grid1, spec_grid, nrow=1, rel_widths = c(1,9), align = "h", axis="tb")
  return(plot_obj)
}

#' Calculate basic statistics summarizing MIMOSA results
#'
#' @import data.table
#' @param node_data Processed table of core MIMOSA results
#' @param threshold_class Whether to use the more or less strict significance threshold for counting metabolite outcomes (default less strict)
#' @return List of counts describing 1) total metabolites analyzed, 2) number consistent, 3) percent consistent, 4) number contrasting, 5) percent contrasting
#' get_summary_stats(node_data)
#' @export
get_summary_stats = function(node_data, threshold_class = "low"){
  totalAnalyzed = node_data[,length(unique(compound))]
  if(threshold_class == "low"){
    totConsistent = node_data[PredictionType2=="Consistent",length(unique(compound))]
    percConsistent = round(node_data[,length(unique(compound[PredictionType2=="Consistent"]))/length(unique(compound))],3)*100
    totContrasting = node_data[PredictionType2=="Contrasting",length(unique(compound))]
    percContrasting = round(node_data[,length(unique(compound[PredictionType2=="Contrasting"]))/length(unique(compound))],3)*100
  } else {
    totConsistent = node_data[PredictionType=="Consistent",length(unique(compound))]
    percConsistent = round(node_data[,length(unique(compound[PredictionType=="Consistent"]))/length(unique(compound))],3)*100
    totContrasting = node_data[PredictionType=="Contrasting",length(unique(compound))]
    percContrasting = round(node_data[,length(unique(compound[PredictionType=="Contrasting"]))/length(unique(compound))],3)*100
  }
  return(c(totalAnalyzed = totalAnalyzed, totConsistent = totConsistent, percConsistent = percConsistent, totContrasting = totContrasting, percContrasting = percContrasting))
}


#' Convert edge list and node list to class network object for plotting with ggnetwork
#'
#' @import data.table
#' @import network
#' @param gene_contribs output of gene_contributions
#' @param node_data output of run_all_metabolites
#' @param node_attrs vector of node attributes to include in network
#' @param edge_attrs vector of edge attributes to include in network
#' @return network object
#' @examples
#' make_contrib_network(gene_contribs_all, node_data, c("Dataset", "Correlation", "Metabolite", "PredictionType2", "SuperPath"), c("KO", "Cor", "Dataset", "stoichReac", "stoichProd"))
#' @export
make_contrib_network = function(gene_contribs, node_data, node_attrs, edge_attrs){
  vertices = unique(c(gene_contribs[,Reac], gene_contribs[,Prod]))
  gene_contribs[,ReacNum:=factor(Reac, levels = vertices)] #Separate by dataset
  gene_contribs[,ProdNum:=factor(Prod, levels = vertices)]
  edge_data = unique(gene_contribs[,list(as.numeric(ReacNum),as.numeric(ProdNum))])#, KO, Cor, stoichReac, stoichProd, Dataset)]) #,compound, Dataset)
  node_data[,nodeNum:=factor(compound, levels = vertices)]
  net1 = network(edge_data, directed = T)

  #Set relevant node attributes
  for(attr in node_attrs){
    if(node_data[,is.factor(get(attr))]){
      set.vertex.attribute(net1, attrname = attr, value = node_data[,as.character(get(attr))], v = node_data[,as.numeric(nodeNum)])
    }else{
      set.vertex.attribute(net1, attrname = attr, value = node_data[,get(attr)], v = node_data[,as.numeric(nodeNum)])
    }
  }
  # all_vertex_ids = unlist(lapply(net1$val, function(x){ return(x$vertex.names)}))
  # all_edge_ids = unique(unlist(sapply(all_vertex_ids, get.edgeIDs, x=net1)))
  for(edge_attr in edge_attrs){
    for(j in 1:nrow(gene_contribs)){
      edgeID = unlist(get.dyads.eids(net1, gene_contribs[j,as.numeric(ReacNum)], gene_contribs[j,as.numeric(ProdNum)]))
      if(!edge_attr %in% names(net1$mel[[edgeID]]$atl)){
        set.edge.attribute(net1, attrname = edge_attr, e=edgeID, value = gene_contribs[j,get(edge_attr)])
      } else{
        old_edge_attr = get.edge.attribute(net1, edge_attr, unlist = F)[[edgeID]]
        if(is.character(old_edge_attr)){
          new_edge_attr = paste0(old_edge_attr, " ", gene_contribs[j,get(edge_attr)])
        } else {
          new_edge_attr = paste0(round(old_edge_attr, 3), " ", round(gene_contribs[j,get(edge_attr)]))
        }
        set.edge.value(net1, attrname = edge_attr, e=edgeID, value = new_edge_attr)
      }
    }
  }
  return(net1)
}

#' Use ggnetwork to make a network visualization of metabolite outcomes and primary gene/reaction contributors
#'
#' @import data.table
#' @import network
#' @import ggplot2
#' @import ggnetwork
#' @param gene_contribs output of gene_contributions
#' @param node_data output of run_all_metabolites
#' @param node_attrs vector of node attributes to include in network
#' @param edge_attrs vector of edge attributes to include in network
#' @return network object
#' @examples
#' plot_contrib_net(gene_contribs_all, node_data, c("Dataset", "Correlation", "Metabolite", "PredictionType2", "SuperPath"), c("KO", "Cor", "Dataset", "stoichReac", "stoichProd"))
#' @export
plot_contrib_net = function(net_obj, col_attr, node_size_attr, edge_size_attr, node_lab, edge_lab){
  net_plot = ggplot(ggnetwork(net_obj,layout = "fruchtermanreingold", cell.jitter = 0.75, niter=1000, repulse.rad = 16000), aes(x=x, y = y, xend=xend, yend=yend)) + geom_edges(arrow = arrow(length = unit(6, "pt"), type = "closed")) + geom_nodes(aes_string(col = col_attr, size = node_size_attr)) + theme_blank() + geom_nodetext_repel(aes_string(label=node_lab), size=3) + geom_edgetext_repel(aes_string(label = edge_lab), size=2, label.padding = unit(0.1, "lines"))
  return(net_plot)
}


#Square plot
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
