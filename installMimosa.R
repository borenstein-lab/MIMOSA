## Install dependencies from CRAN and Bioconductor
if(!requireNamespace("devtools", quietly = T)){
 install.packages("devtools", repos = "https://cloud.r-project.org/")
}
cran_packages = c("devtools", "data.table", "Rcpp", "getopt", "RColorBrewer", "rmarkdown", "BiocManager", "vegan", "permute", "network", "ggnetwork", "cowplot", "readr", "testthat", "pandoc")
for (package in cran_packages){
 if(!requireNamespace(package, quietly = T)){
   if(package=="vegan"){ #Vegan not currently working with devtools for unknown reason
     install.packages("vegan", repos = "https://cloud.r-project.org/")
   }
   devtools::install(package)
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
