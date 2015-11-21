# Jinliang Yang
# Oct. 29, 2015

# install packages
#devtools::install_github("hadley/devtools")
#library(devtools) 
#install_github("vsbuffalo/tasselr") 
#install_github("vsbuffalo/ProgenyArray")


library("tasselr")
library("imputeR")

# Note: at least 100G memory was needed to load the hdf5 file
# load h5file
land <- initTasselHDF5("largedata/maize_landrace.h5", version="5")
land <- loadBiallelicGenotypes(land, verbose = TRUE)

# reformat to imputeR object
ob2 <- imputeRob(h5=land, missingcode=3)
save(file="largedata/bode_data.RData", list="ob2")
