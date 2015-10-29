# jinliang Yang

library(devtools) 
install_github("vsbuffalo/tasselr") 
install_github("vsbuffalo/ProgenyArray")



library(parallel)
options(mc.cores=NULL)
# you need to specify the location where the packages were installed. 
load_all("~/bin/tasselr")
load_all("~/bin/ProgenyArray")
load_all("~/Documents/Github/imputeR")

# Note: at least 64G memory was needed to load the hdf5 file
ob <- loading_hdf5(h5file="largedata/teo.h5")