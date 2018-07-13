# MIMOSA
Model-based Integration of Metabolite Observations and Species Abundances

A pipeline for joint metabolic model-based analysis of metabolomics measurements and taxonomic composition from microbial communities. Includes an R package that implements the core integrative analysis functions as well as example pre-processing and analysis scripts. MIMOSA is under active development, but version 1.0.0 is a general implementation of the methods described and used in (Noecker et al, mSystems 2016). See http://elbo.gs.washington.edu/download.html for a previous version of the code initially used to produce the results in that work.

## Contents

- **mimosa**: An R package to generate community-specific metabolic network model, use the model to make metabolite predictions and identify consistent and contrasting metabolite variation, 
identify potential taxonomic and gene contributors to metabolite variation, and shuffle the community metabolic network for comparison of results with a null model. 
The core analysis uses the file `reaction_mapformula.lst` from the KEGG database (Kanehisa and Goto, 2000) to link genes to reactions.

You can install **mimosa** using the `devtools` package in R:
```R
devtools::install_github("borenstein-lab/MIMOSA/mimosa")
```

If the above command produces an error related to installing dependencies, it most likely indicates that you need to install the **qvalue** or **KEGGREST** from [Bioconductor](https://bioconductor.org/packages/release/bioc/html/qvalue.html), as these are not available from CRAN.

- **runMimosa.R**: Script used to run main MIMOSA analyses from the command line (see "Using MIMOSA").

- **run_all.sh**: example usage of `runMimosa.R` to regenerate results from one of the datasets described in Noecker et al, 2016. You can download these files from [here](http://elbo.gs.washington.edu/download.html). This script uses the Python modules [`picrust`](http://picrust.github.io/picrust/)(Langille et al, 2013) and [`MUSiCC`](http://elbo.gs.washington.edu/software_musicc.html)(Manor and Borenstein, 2015) to process and prepare the datasets. 

- **summarizeMIMOSAresults.Rmd** and **processMimosaOut.R**: example scripts for summarizing and plotting the output of MIMOSA.

## Using MIMOSA

You can run a full MIMOSA analysis by running the script runMimosa.R from the command line with additional arguments and the output will be saved to multiple files. You can also run any of the individual steps in R. 

You can see example data file formatting in the `tests/testthat` directory. Gene and metabolite abundance files should be formatted like the example files, and sample IDs need to be consistent between the two files.

The analysis consists of the following main steps:

1. Construct a community-specific metabolic network model linking genes, reactions, and metabolites.

2. Use the above model to calculate community metabolic potential scores for each sample and metabolite and compare those using Mantel tests with metabolite concentrations.

3. Identify using a heuristic approach which genes and reactions were important contributors to variation in community metabolic potential scores.

4. Identify taxa that are consistent with contributing to metabolic potential scores.

5. Optionally, perform a permutation test of the metabolic network to quantify how much the network model strengthens associations between the two data types.


### Required arguments for runMimosa.R

- **-g,--genefile**: file path to gene (KO) abundances (total across all species), preferably normalized with MUSiCC. Example: `tests/testthat/test_genes.txt`
- **-m,--metfile**: file path to metabolite abundances (with KEGG compound IDs), preferably unscaled peak abundances or absolute concentrations. Example: `tests/testthat/test_mets.txt`
- **-o,--contribs_file**: file path to metagenome contributions (species-specific gene abundances), formatted like the output of the PICRUSt `metagenome_contributions.py` script. Example: `tests/testthat/test_contributions.txt`
- **-e,--mapformula_file**: File path to KEGG reaction_mapformula.lst file (included in FTP download of KEGG database), containing information on reactions annotated in KEGG pathways. Example: `tests/testthat/test_mapformula.txt`
- **-p,--file_prefix**: file path prefix where output files will be written


### Optional arguments for runMimosa.R

The following two files are required for full reaction information and stoichiometry, if you are using a downloaded version of KEGG. If they are not provided, MIMOSA will try to use the KEGGREST package to get the same information from the KEGG API.
- **-r,--ko_rxn_file** File path to links between KO gene families and reactions (found at genes/ko/ko_reaction.list in the KEGG FTP download)
- **-x,--rxn_annots_file** File path to full KEGG Reaction annotations (found at ligand/reaction/reaction in the KEGG FTP download)

Other optional arguments:

- **-b,--keggFile**: file path to pre-computed generic community network template (output of the function `generate_network_template_kegg`), if you have already generated that information
- **-u,--num_permute**: Number of permutations to perform for the Mantel test comparison (default: 20000)
- **-c,--cor_method**: Whether to use Spearman or Pearson correlations for the Mantel test (default and recommended: Spearman)
- **-f,--degree_filter**: Connectivity threshold filter. Metabolites connected in the community network model to this number of KOs or higher are considered currency metabolites and are filtered from the analysis (default: 30).
- **-z,--nonzero_filt**: Metabolites found to have a nonzero concentration or nonzero metabolic potential scores in this number of samples or fewer are filtered from the analysis (default: 3)
- **-t,--taxonomy_file**: File path to taxonomic information for each OTU. Include if you would like to evaluate contributions at the genus level rather than the OTU level. The file must include a column named "OTU" that matches the OTUs in the contributions file. It can also include either a column named "Genus" or a second column that contains a full 7-level taxonomy in the Greengenes format (e.g. "k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Prevotellaceae; g__Prevotella; s__")
- **-d, --metadata_file**: File path to sample metadata file (used when generating final summary of MIMOSA results to compare cases-enriched vs control-enriched metabolites). To include this analysis, this final must contain a table with a column named "Sample" containing sample IDs and a column named with the metadata variable (below) specifying cases and controls as 1s and 0s.
- **-v, --metadata_var**: Variable column name (from metadata file) to use as binary indicator to calculate metabolite enrichment when generating final summary of MIMOSA results

### Example usage of runMimosa.R

**Option 1: Using downloaded KEGG files to build community network**
```bash
Rscript runMimosa.R --genefile="Dataset2_picrust_musicc.txt" --metfile="Dataset2_mets.txt" --contribs_file="Dataset2_metagenome_contributions.txt" --file_prefix="Dataset2_bv" --mapformula_file="keggPath/reaction_mapformula.lst" --ko_rxn_file="keggPath/ko_reaction.list" --rxn_annots_file="keggPath/reaction" --nonzero_filt=4 
```

**Option 2: Using the KEGGREST package to build the community metabolic network template (can be slow depending on network connection)**
```bash
Rscript runMimosa.R --genefile="Dataset2_picrust_musicc.txt" -m "Dataset2_mets.txt" --contribs_file="Dataset2_metagenome_contributions.txt" --file_prefix="Dataset2_bv" --mapformula_file="keggPath/reaction_mapformula.lst" --nonzero_filt=4 
```

**You can test that the package is installed and working using the test data.** 

```bash
Rscript runMimosa.R --genefile="mimosa/tests/testthat/test_genes.txt" --metfile="mimosa/tests/testthat/test_mets.txt" --contribs_file="mimosa/tests/testthat/test_contributions.txt" --mapformula_file="mimosa/tests/testthat/test_mapformula.txt" --file_prefix="test1" --ko_rxn_file="mimosa/tests/testthat/test_ko_reaction.txt" --rxn_annots_file="mimosa/tests/testthat/test_reaction.txt" --metadata_file="mimosa/tests/testthat/test_metadata.txt" --metadata_var="BV"
```


### Output of runMimosa.R

The main MIMOSA results will be summarized in a series of plots and tables in the file **summarizeMIMOSAresults.html**. Other results files with full details are specified below.

Core comparison analysis:
- **file_prefix_nodes.txt**: A table of all successfully analyzed metabolites, their Mantel correlation between metabolic potential scores and concentrations, and the associated p-values and q-values (FDR corrected).
- **file_prefix_edges.txt**: The full community network model for that dataset, with gene, reaction, product, and stoichiometry information.
- **file_prefix_signifPos.txt/file_prefix_signifNeg.txt**: subsets of the "nodes" file that were significant at an FDR of 0.1.
- **file_prefix_cmpAll.txt**: Table of metabolic potential scores across all samples and metabolites.
- **file_prefix_out.rda**: all of output of the run_all_metabolites for the core comparison analysis that can be loaded into an R session.
- **full_network_template_all_info.txt**: Table of the full KEGG community network model template combining information on gene families, reactions, and stoichiometry.

Potential taxonomic contributor analysis:
- **file_prefix_specContrib.txt**: Table of results of analysis of potential taxonomic contributors for each metabolite (correlations of single-species scores with full community score).
- **file_prefix__allKOAbundsByOTU.rda, file_prefix_allCMPsAloneByOTU.rda**: Intermediate Rdata files with gene abundances and metabolic potential scores for each species alone.


Potential gene/species contributor analysis:
- **file_prefix_geneContribAnalysis.txt**: Table of results of contribution analysis of most relevant potential gene and reaction contributors for each metabolite (correlations indicate how much removing the reactions associated with that gene affected metabolite potential scores).
- **file_prefix_geneContribCompoundSummary.txt**: Table summarizing gene contribution results by compound - number of relevant synthesizing and degrading genes for each metabolite and whether it was primarily predicted based on synthesis, degradation, or a combination.
- **file_prefix_geneContribs.rda**: full output of gene contribution analysis.


## Citation

Noecker, C., Eng, A., Srinivasan, S., Theriot, C.M., Young, V.B., Jansson, J.K., Fredricks, D.N., and Borenstein, E. (2016). Metabolic Model-Based Integration of Microbiome Taxonomic and Metabolomic Profiles Elucidates Mechanistic Links between Ecological and Metabolic Variation. mSystems 1, e00013–15, doi:[10.1128/mSystems.00013-15](http://msystems.asm.org/content/1/1/e00013-15).
