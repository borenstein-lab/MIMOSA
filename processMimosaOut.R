##Make a table of consistent metabolites and relevant info:
options(stringsAsFactors = F)
library(data.table)
library(mimosa)
args = commandArgs(trailingOnly = T)
prefix = args[1]
contribs_file = args[2]

consistent_node_data = get_all_consistent_metabolite_info(prefix, contribs_file)

write.table(consistent_node_data, file = paste0(prefix, "_consistentPos_withContributors.txt"), quote=F, row.names = F, sep = "\t")

#Now every row in consistent_node_data should represent genes and species that were identified as potential important contributors for a consistent metabolite so you can go through this and see what looks interesting, count how many different metabolites each species was important for, etc.
