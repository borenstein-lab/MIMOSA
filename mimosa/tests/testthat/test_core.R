library(data.table)
options(stringsAsFactors = F)
context("Core MIMOSA tests")

test_gene_file = "test_genes.txt"
test_met_file = "test_mets.txt"
test_ko_rxn_file = "test_ko_reaction.txt"
test_rxns_file = "test_reaction.txt"
test_mapformula_file = "test_mapformula.txt"
test_contrib_file = "test_contributions.txt"

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

all_kegg = get_kegg_reaction_info(test_ko_rxn_file, test_rxns_file, save_out = F, kolist = datasets[[1]][,KO])
rxn_table = generate_network_template_kegg(test_mapformula_file, all_kegg, write_out = F)

test_that("Network data looks normal", {
  #expect_gt(nrow(get_kegg_reaction_info("KEGGREST", kolist = datasets[[1]][,KO])),1)
  expect_equal(names(all_kegg), c("KOs", "Reactions", "Reaction_info", "kos_to_rxns"))
  expect_lte(rxn_table[,length(unique(KO))], datasets[[1]][,length(KO)])
  expect_gt(rxn_table[,length(unique(KO))], 0)
  expect_equal(names(rxn_table), c("Rxn", "KO", "Reac", "Prod", "Path", "ReacProd", "stoichReac", "stoichProd"))
})

ko_net = generate_genomic_network(datasets[[1]][,KO], keggSource = "KeggTemplate", rxn_table = rxn_table)
ko_net2 = generate_genomic_network(datasets[[1]][,KO], keggSource = "KeggTemplate", rxn_table = rxn_table, normalize = F)

test_that("Networks generated successfully", {
  expect_length(ko_net, 3)
  expect_equal(names(ko_net[[3]]), c("KO", "Reac", "Prod", "stoichReac", "stoichProd"))
  expect_equal(sort(unique(names(ko_net[[1]]))), sort(unique(ko_net[[3]][,KO])))
  expect_equal(sort(unique(row.names(ko_net[[1]]))), sort(unique(c(ko_net[[3]][,Prod], ko_net[[3]][,Reac]))))
  expect_length(ko_net2, 3)
  expect_equal(names(ko_net2[[3]]), c("KO", "Reac", "Prod", "stoichReac", "stoichProd"))
  expect_equal(sort(unique(names(ko_net2[[1]]))), sort(unique(ko_net2[[3]][,KO])))
  expect_equal(sort(unique(row.names(ko_net2[[1]]))), sort(unique(c(ko_net2[[3]][,Prod], ko_net2[[3]][,Reac]))))
})

cmp_scores = get_cmp_scores(ko_net[[1]], datasets[[1]])
cmp_scores2 = get_cmp_scores(ko_net2[[1]], datasets[[1]])

test_that("CMP scores can be calculated from network", {
  expect_gt(nrow(cmp_scores), 0)
  expect_equal(ncol(cmp_scores), ncol(datasets[[1]]))
  expect_gt(nrow(cmp_scores2), 0)
  expect_equal(ncol(cmp_scores2), ncol(datasets[[1]]))
})

shared_mets = intersect(datasets[[2]][,KEGG], cmp_scores[,compound])
met_mat = make_pairwise_met_matrix(metabolite = shared_mets[1], met_mat = datasets[[2]])
met_mat2 = make_pairwise_met_matrix(shared_mets[1], cmp_scores)
met_mat_a = make_pairwise_met_matrix(shared_mets[4], datasets[[2]])
met_mat_2a = make_pairwise_met_matrix(shared_mets[4], cmp_scores)

test_that("Cmp scores can be compared with mets", {
  expect_gt(nrow(met_mat), 0)
  expect_equal(nrow(met_mat), ncol(met_mat))
  expect_equal(nrow(met_mat2), ncol(met_mat2))
  expect_equal(nrow(met_mat), nrow(met_mat2))
  expect_error(mantel_2sided(met_mat, met_mat2, permutations = 500, direction = "pos", method = "spearman"))
  expect_silent(mantel_2sided(met_mat_a, met_mat_2a, permutations = 500, direction = "pos", method = "spearman"))
})

run_all_metabolites(datasets[[1]], datasets[[2]], file_prefix = "test", id_met = F,
                    net_method = "KeggTemplate", net_file = net_file, rxn_table_source = rxn_table,
                    correction = "fdr", degree_filter = 30, cor_method = "spearman", nperm = 200, nonzero_filter = 4)
load("test_out.rda")

test_that("Run all metabolites function works", {
  #expect_equal(length(all_comparisons), length(shared_mets))
  expect_equal(nrow(node_data), length(all_comparisons))
})

spec_contribs = get_spec_contribs(test_contrib_file, data_dir = getwd(), results_file = "test_out.rda", out_prefix = "test", otu_id = "all", valueVar = "singleMusicc",
                                  sum_to_genus = T, write_out = T, taxonomy_file = "test_taxonomy.txt")
spec_contribs2 = get_spec_contribs(test_contrib_file, data_dir = getwd(), results_file = "test_out.rda", out_prefix = "test", otu_id = "all", valueVar = "RelAbundSample", sum_to_genus = T, write_out = T, taxonomy_file = "test_taxonomy.txt")


test_that("Species contributions work", {
  expect_gt(nrow(spec_contribs), 0)
  expect_gt(nrow(spec_contribs[!is.na(Cor)]), 0)
  expect_gt(nrow(spec_contribs2), 0)
  expect_gt(nrow(spec_contribs2[!is.na(Cor)]), 0)
})
