## Install dependencies from CRAN and Bioconductor

cran_packages = c("devtools", "data.table", "Rcpp", "getopt", "RColorBrewer", "rmarkdown", "BiocManager")
for (package in cran_packages){
 if(!requireNamespace(package, quietly = T)){
   install.packages(package)
 }
}

bioc_packages = c("qvalue", "KEGGREST")
for (package in bioc_packages){
 if(!requireNamespace(package, quietly = T)){
  BiocManager::install(package)
 }
}

## Install MIMOSA
install.packages("mimosa/", repos=NULL, type = "source")
