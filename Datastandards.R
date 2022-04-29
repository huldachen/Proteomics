#mzR:provides a unified interface to various mass spectrometry open formats
#Load the required packages
library("mzR") #Software package
library("msdata") #Data package
#Extract the releavant example file from the local 'msdata' installation
filepath <- system.file("microtofq", package = "msdata")
file <- list.files(filepath, pattern="MM14.mzML",
                   full.names=TRUE, recursive = TRUE)
#Create a commection to the mzML file
mz <- openMSfile(file)
#Demonstrate data access
basename(fileName(mz))
runInfo(mz)
instrumentInfo(mz)
#It's better to close the connection explicitely
close(mz)

## Identification file handle.
file <- system.file("mzid", "Tandem.mzid.gz", package="msdata")
mzid <- openIDfile(file)
softwareInfo(mzid)
enzymes(mzid)
names(psms(mzid))
head(psms(mzid))[, 1:13]


# MSnbase:base functions and classes for MS-based proteomics
library("MSnbase")

#Simple test included in the package
mzXML <- dir(system.file(package="MSnbase",dir="extdata"),
             full.name=TRUE,
             pattern="mzXML$")
basename(mzXML)

#Reads the raw data into and MSnExp instance
raw <- readMSData(mzXML, verbose = FALSE, centroided = TRUE)
raw

#Extract a single spectrum
raw[[3]]
plot(raw, full = TRUE)
plot(raw[[3]], full = TRUE, reporters = iTRAQ4)