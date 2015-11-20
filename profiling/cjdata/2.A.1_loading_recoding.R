# Jinliang Yang
# Oct. 29, 2015

library("tasselr")
library("imputeR")

# Note: at least 100G memory was needed to load the hdf5 file
# load h5file
teod <- initTasselHDF5("largedata/teo.h5", version="5")
teod <- loadBiallelicGenotypes(teod, verbose = TRUE)

# reformat to imputeR object
ob <- imputeRob(h5=teod, missingcode=3)
save(file="largedata/teo.RData", list="ob")

#calculating missing rates ... done.
#calculating minor allele frq (MAF) ... done.
#calculating reference allele frq (RAF) ... done.
#preparing Geno4imputeR object for [4875] plants and [598043] SNPs ... 