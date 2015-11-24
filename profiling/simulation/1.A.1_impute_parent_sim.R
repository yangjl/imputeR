# simulaltion 
set.seed(1234)
GBS.array <- sim.array(size.array=50, numloci=100, hom.error = 0.02, het.error = 0.8,
                       rec = 0.25, selfing = 0.5, imiss = 0.5, misscode = 3)

We now impute parent genotypes using `impute_parent` and extract results using `parentgeno`. In the resulting table `res`, the first three columns are the probabilities of genotype `0, 1, 2`. The 4th column is the odd ratio of the highest divided by the 2nd highest probability. `gmax` parent's genotype with the highest probability. `gor` parent's genotype with the highest probability and `OR` bigger than the specified threshold.

```{r, eval=TRUE}
#The stuff:
inferred_geno_likes <- impute_parent(GBS.array, hom.error=0.02, het.error=0.8, imiss=0.5, p = GBS.array@freq)
res <- parentgeno(inferred_geno_likes, oddratio=0.6931472, returnall=TRUE)
res$true_parent <- GBS.array@true_parents[[50]]$hap1 + GBS.array@true_parents[[50]]$hap2
#error rates
nrow(subset(res, gmax != true_parent ))/nrow(res)
nrow(subset(res, gor != true_parent & gor !=3 ))/nrow(res)