# Jinliang Yang
# Oct. 29, 2015

library("tasselr")
library("imputeR")

# Note: at least 100G memory was needed to load the hdf5 file
# load h5file
teo <- initTasselHDF5("largedata/teo.h5", version="5")
teo <- loadBiallelicGenotypes(teo, verbose = TRUE)

# reformat to imputeR object
ob1 <- imputeRob(h5=teo, missingcode=3)
save(file="largedata/teo.RData", list="ob1")

#calculating missing rates ... done.
#calculating minor allele frq (MAF) ... done.
#calculating reference allele frq (RAF) ... done.
#preparing Geno4imputeR object for [4875] plants and [598043] SNPs ... 