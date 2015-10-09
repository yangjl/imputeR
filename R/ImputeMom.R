#'
#' \code{Impute mom's genotype using mom and kids' GBS data. } 
#'
#' We have obs. mom and obs. (selfed) kids. We want to know the most likelihood of the mom's genotype. 
#' This function is to impute mom's genotype from a progeny array of k kids at a single locus.
#' inferred_mom=1 -> 00, 2->01, 3->11
#'
#' @param obs_mom A vector of mom's GBS data. Coded with 0, 1 and 2, the copy of reference alleles. 
#' Missing data should be coded with 3.
#' @param obs_kids A list of vectors of Kid's GBS data. Coded with 0, 1 and 2, the copy of reference alleles. 
#' Missing data should be coded with 3.
#' @param  p p=sample(sfs,numloci) #freqs of all loci sampled from site frequency spectrum (SFS)
#' @return If all inputs are integer and logical, then the output
#'   will be an integer. If integer overflow
#'   \url{http://en.wikipedia.org/wiki/Integer_overflow} occurs, the output
#'   will be NA with a warning. Otherwise it will be a length-one numeric or
#'   complex vector.
#'
#'   Zero-length vectors have sum 0 by definition. See
#'   \url{http://en.wikipedia.org/wiki/Empty_sum} for more details.
#' @examples
#' sum(1:10)
#' sum(1:5, 6:10)
#' sum(F, F, F, T, T)
#'
#' sum(.Machine$integer.max, 1L)
#' sum(.Machine$integer.max, 1)
#'
#' \dontrun{
#' sum("a")
#' }

impute_mom <- function(obs_mom, obs_kids, locus,  p){
    ### locus:
    ### p: p=sample(sfs,numloci) #freqs of all loci sampled from sfs
    mom_probs <- as.numeric()
    for(inferred_mom in 1:3){
        #P(G'|G)
        pgg <- gen_error_mat[inferred_mom, obs_mom[locus]+1] #+1 because obs_mom is 0,1, or 2
        #P(G)
        pg <- hw_probs(p[locus])[inferred_mom]
        #P(kids|G) sum of logs instead of product
        pkg <- sum(sapply(1:length(obs_kids), function(z){
            ### take care of missing data
            #if(progeny[[z]][[2]][locus] >= 0 & progeny[[z]][[2]][locus] <=2){
            log(sum(probs[[inferred_mom]][, obs_kids[[z]][[2]][locus]+1]))
            #}
        } ))
        mom_probs[inferred_mom] <- pkg+log(pgg)+log(pg)
    }
    return(which(mom_probs==max(mom_probs))-1)
}
############################################################################

## Functions
############################################################################
# Return HW probs
hw_probs <- function(x){ return(c(x^2,2*x*(1-x),(1-x)^2))}
############################################################################
# row 1 is true_gen 00, row2 is true_gen 01, row 3 is true_gen 11
# cols are obs. genotype (00,01,11)

get_error_mat <- function(hom.error, het.error){
    gen_error_mat <- matrix(c(1-hom.error,hom.error/2,hom.error/2,het.error/2,
                            1-het.error,het.error/2,hom.error/2,hom.error/2,1-hom.error),
                            byrow=T,nrow=3,ncol=3)
    
    probs <- vector("list",3)
    ### 1: 00x01; 2: 01x01; 3: 11x01
    probs[[1]] <- gen_error_mat*matrix(c(1/2, 1/2, 0), nrow = 3,byrow=F,ncol=3)
    probs[[2]] <- gen_error_mat*matrix(c(1/4, 1/2, 1/4), nrow = 3,byrow=F,ncol=3)
    probs[[3]] <- gen_error_mat*matrix(c(0, 1/2, 1/2), nrow = 3,byrow=F,ncol=3)
    
    gen_error_mat <- cbind(gen_error_mat, 1)
    probs[[1]] <- cbind(probs[[1]], 1)
    probs[[2]] <- cbind(probs[[2]], 1)
    probs[[3]] <- cbind(probs[[3]], 1)
    
    return(list(gen_error_mat, probs))
}


