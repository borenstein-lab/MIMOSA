library(data.table)
options(stringsAsFactors = F)
context("Core MIMOSA tests")

test_gene_file = "test_genes.txt"
test_met_file = "test_mets.txt"
test_ko_rxn_file = "test_ko_reaction.txt"
test_rxns_file = "test_reaction.txt"
test_mapformula_file = "test_mapformula.txt"

datasets = read_files(test_gene_file, test_met_file)

test_that("File reading results in no NAs and consistent naming", {
  expect_equal(ncol(datasets[[1]]), ncol(datasets[[2]]))
  expect_equal(length(intersect(names(datasets[[1]]), names(datasets[[2]]))), ncol(datasets[[1]])-1)
  expect_gt(nrow(datasets[[1]]), 0)
  expect_gt(nrow(datasets[[2]]), 0)
  expect_gt(ncol(datasets[[1]]), 1)
  expect_false(any(is.na(datasets[[1]])))
  expect_false(any(is.na(datasets[[2]])))
})

test_that("Network can be generated using all potential source data", {
  #expect_gt(nrow(get_kegg_reaction_info("KEGGREST", kolist = datasets[[1]][,KO])),1)
  all_kegg = get_kegg_reaction_info(test_ko_rxn_file, test_rxns_file, save_out = F, kolist = datasets[[1]][,KO])
  expect_equal(names(all_kegg), c("KOs", "Reactions", "Reaction_info", "kos_to_rxns"))
  rxn_table = generate_network_template_kegg(test_mapformula_file, all_kegg, save_out = F)
  expect_lte(rxn_table[,length(unique(KO))], datasets[[1]][,length(KO)])
  expect_gt(rxn_table[,length(unique(KO))], 0)
  expect_equal(names(rxn_table), c("Rxn", "KO", "Reac", "Prod", "Path", "ReacProd", "stoichReac", "stoichProd"))
  ko_net = generate_genomic_network(datasets[[1]][,KO], keggSource = "KeggTemplate", rxn_table = rxn_table)
  expect_length(ko_net, 3)
  expect_equal(names(ko_net[[3]]), c("KO", "Reac", "Prod", "stoichReac", "stoichProd"))
  expect_equal(sort(unique(names(ko_net[[1]]))), sort(unique(ko_net[[3]][,KO])))
  expect_equal(sort(unique(row.names(ko_net[[1]]))), sort(unique(c(ko_net[[3]][,Prod], ko_net[[3]][,Reac]))))
})

