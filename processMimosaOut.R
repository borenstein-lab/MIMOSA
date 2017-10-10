##Make a table of consistent metabolites and relevant info (in R):
options(stringsAsFactors = F)
library(data.table)
library(mimosa)
args = commandArgs(trailingOnly = T)
prefix = args[1]
contribs_file = args[2]

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

write.table(consistent_node_data, file = paste0(prefix, "_consistentPos_withContributors.txt"), quote=F, row.names = F, sep = "\t")
#Now every row in consistent_node_data should represent genes and species that were identified as potential important contributors for a consistent metabolite so you can go through this and see what looks interesting, count how many different metabolites each species was important for, etc.
