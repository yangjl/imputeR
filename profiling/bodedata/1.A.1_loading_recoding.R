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
load_all("~/Documents/Github/imputeR")

# Note: at least 100G memory was needed to load the hdf5 file
# load h5file
land <- initTasselHDF5("largedata/maize_landrace.h5", version="5")
land <- loadBiallelicGenotypes(land, verbose = TRUE)

# reformat to imputeR object
source("R/my-TasselHDF5-method.R")
ob <- imputeRob(land)
save(file="largedata/bode_data.RData", list="ob")
