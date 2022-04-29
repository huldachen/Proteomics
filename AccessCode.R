# gets the vignette source
rfile <- system.file("doc/RforProteomics.R",
                     package = "RforProteomics")
rfile

#Prepare the working environment
library("RColorBrewer") # Color palettes
library("ggplot2")  # Convenient and nice plotting
library("reshape2") # Flexibly reshape data