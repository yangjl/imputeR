### Jinliang Yang
### 12/5th, 2015

library(imputeR)

# to make the random events repeatable
set.seed(12345)
# simulate a GBS.array object
GBS.array <- sim.array(size.array=30, numloci=1000, hom.error = 0.02, het.error = 0.8,
                       rec = 0.25, selfing = 0, imiss = 0.5, misscode = 3)
# get perfect parent genotype
GBS.array <- get_true_GBS(GBS.array, phased.parents = TRUE)
# get probability matrices
probs <- error_mx2(hom.error=0.02, het.error=0.8)
# phasing   
phase <- phase_parent(GBS.array, win_length=10, join_length=10, 
                      self_cutoff = 100, verbose=TRUE)

# compute error rate
out <- phase_error_rate(GBS.array, phase)

outfile <- paste0("largedata/sim2/size", SIZE, "_oc_10kloci.csv")
write.table(out, outfile, sep=",", row.names=FALSE, quote=FALSE)



