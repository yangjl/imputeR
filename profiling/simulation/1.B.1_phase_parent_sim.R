### Jinliang Yang
### 12/5th, 2015

library(devtools)
load_all()

# to make the random events repeatable
set.seed(123457)
# simulate a GBS.array object
GBS.array <- sim.array(size.array=40, numloci=100, hom.error = 0.02, het.error = 0.8,
                       rec = 0.25, selfing = 1, imiss = 0.5, misscode = 3)
# get perfect parent genotype
GBS.array <- get_true_GBS(GBS.array, phased.parents = TRUE, parents.mr=0.1)
# get probability matrices
probs <- error_mx2(major.error=0.02, het.error=0.8, minor.error=0.02)
# phasing   
phase <- phase_parent(GBS.array, win_length=10, join_length=10, verbose=TRUE, OR=log(3))

# compute error rate
out <- phase_error_rate(GBS.array, phase)

outfile <- paste0("largedata/sim2/size", SIZE, "_oc_10kloci.csv")
write.table(out, outfile, sep=",", row.names=FALSE, quote=FALSE)



