# Jinliang Yang
# Oct. 29, 2015

# install packages
#devtools::install_github("hadley/devtools")
#library(devtools) 
#install_github("vsbuffalo/tasselr") 
#install_github("vsbuffalo/ProgenyArray")


# load packages
library(parallel)
library(devtools)
options(mc.cores=NULL)
# you need to specify the location where the packages were installed. 
load_all("~/bin/tasselr")
load_all("~/bin/ProgenyArray")
load_all("~/Documents/Github/imputeR")

# Note: at least 100G memory was needed to load the hdf5 file
# load h5file
teod <- initTasselHDF5("largedata/teo.h5", version="5")
teod <- loadBiallelicGenotypes(teod, verbose = TRUE)
save(file="largedata/teod.RData", list="teod")

# reformat to imputeR object
ob <- imputeRob(teod)
save(file="largedata/teo.RData", list="ob")
