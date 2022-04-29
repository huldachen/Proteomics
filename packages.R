#List all the vignettes in the RforProteomics
vignette(package = "RforProteomics")
vignette("RforProteomics", package = "RforProteomics")

#First time need to install Bioconductor packages
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

#Else, load package
library("BiocManager")
BiocManager::install("RforProteomics")

# Some packages used in the document depend on external libraries that need to be installed prior to the R:
