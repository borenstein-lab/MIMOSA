## Install dependencies from CRAN and Bioconductor
if(!requireNamespace("devtools", quietly = T)){
 install.packages("devtools", repos = "https://cloud.r-project.org/")
}
cran_packages = c("devtools", "data.table", "Rcpp", "getopt", "RColorBrewer", "rmarkdown", "BiocManager", "vegan", "permute", "network", "ggnetwork", "cowplot", "readr", "testthat")
for (package in cran_packages){
 if(!requireNamespace(package, quietly = T)){
   devtools::install_cran(package, repos = "https://cloud.r-project.org")
 }
}

bioc_packages = c("qvalue", "KEGGREST")
for (package in bioc_packages){
 if(!requireNamespace(package, quietly = T)){
  BiocManager::install(package)
 }
}

## Install MIMOSA
devtools::install_local("mimosa/", force = T, dependencies = F)
