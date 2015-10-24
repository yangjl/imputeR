#'
#' \code{Impute mom's genotype using mom and kids' GBS data. } 
#'
#' We have observed mom and observed (selfed) kids. We want to know the likelihood of the mom's genotypes. 
#' This function is to impute mom's most likely genotype from a progeny array of k kids by 
#' giving a log likelihood threshold.
#'
#' @param obs_mom A vector of mom's GBS data. Coded with 0, 1 and 2, which is the copy of alternative alleles. 
#' Missing data should be coded with 3.
#' @param obs_kids A list of vectors of Kid's GBS data. Coded with 0, 1 and 2, which is the copy of alternative alleles. 
#' Missing data should be coded with 3.
#' @param hom.error Homozygous error rate, default=0.02.
#' @param het.error Heterozygous error rate, default=0.8.
#' @param p vector of allele frequencies at each locus. should be estimated from set of all parents.
#' WARNING: if not supplied, a random neutral SFS is used to generate p at each locus.
#' @param oddratio The cutoff used for determining the log likelihood ratio of the highest and the 2nd highest genotypes. 
#' The oddratio = NULL means to report the most likely genotype. Default value sets to 0.5.
#' @param returnall For function 'genomom', returnall=TRUE will return with all the information; returnall=FALSE will return
#' mom's genotype only.
#' @return return a data.frame with all the log likelihoods. Using function \code{momgeno} to extract mom's genotype. 
#'   
#'   See \url{https://github.com/yangjl/imputeR/blob/master/vignettes/imputeR-vignette.pdf} for more details.
#' @examples
#' obs_mom <- c(0, 0, 0)
#' obs_kids <- list(c(0, 0, 0), c(0, 0, 0), c(0, 0, 0), c(0, 0, 0), c(0, 0, 0), c(0, 0, 0),
#' c(0, 0, 0), c(0, 0, 0), c(0, 0, 0), c(0, 0, 0), c(0, 0, 0), c(0, 0, 0))
#'
#' impute_mom(obs_mom, obs_kids, hom.error=0.02, het.error=0.8)
#' 
#' momgeno(geno, oddratio=0.5, returnall=FALSE)
#' 
#' 
impute_mom <- function(obs_mom, obs_kids, hom.error=0.02, het.error=0.8,p=NULL){
    
    ### need to check genotypes
    numloci <- length(obs_mom)
    
    message(sprintf("###>>> impute mom's genotype using [ %s ] kids in a full sib family ...", length(obs_kids)))
    ### get probability matrices
    gen_error <- gen_error_mat(hom.error, het.error)
    probs <- error_mx(hom.error, het.error)
    
    ### get sfs if not provided?
    if(is.null(p)){
        message(sprintf("###>>> no allele frequencies provided. generating random allele frequencies from a neutral SFS"))
        sfs <- getsfs(x = 1:99/100)
        p <- sample(sfs, numloci,replace=TRUE) #get freqs for all loci
    }
    
    res <- lapply(1:length(obs_mom), function(locus) impute_one_site(locus, gen_error, p, probs))
    geno <- as.data.frame(matrix(unlist(res), ncol=3, byrow=TRUE))
    names(geno) <- c("g0", "g1", "g2")
    
    
    return(geno)
}

#' @rdname impute_mom
momgeno <- function(geno, oddratio=0.5, returnall=TRUE){
    geno$OR <- apply(geno, 1, function(v){
        n <- length(v)
        return(max(v) - sort(v, partial=n-1)[n-1])  
    })
    geno$gmax <- apply(geno[, 1:3], 1, function(v){
        return(which.max(v)-1)  
    })
    geno$gor <- 3
    geno[geno$OR > oddratio, ]$gor <- geno[geno$OR > oddratio, ]$gmax
    
    if(returnall){
        return(geno)
    }else{
        if(is.null(oddratio)){
            return(geno$g)    
        }else{
            return(geno$gor)
        } 
    }
}

#' @rdname impute_mom
impute_one_site <- function(locus, gen_error, p, probs){
    ### locus:
    ### p: p=sample(sfs,numloci) #freqs of all loci sampled from sfs
    mom_probs <- as.numeric()
    for(inferred_mom in 1:3){
        #P(G'mom|G)
        pgg <- gen_error[inferred_mom, obs_mom[locus]+1] #+1 because obs_mom is 0,1, or 2
        #P(G)
        pg <- hw_probs(p[locus])[inferred_mom]
        #P(kids|G) sum of logs instead of product
        pkg <- sum(sapply(1:length(obs_kids), function(z){
            ### take care of missing data
            #if(progeny[[z]][[2]][locus] >= 0 & progeny[[z]][[2]][locus] <=2){
            log(sum(probs[[inferred_mom]][, obs_kids[[z]][locus]+1]))
            #}
        } ))
        mom_probs[inferred_mom] <- pkg+log(pgg)+log(pg)
    }
    #return(which(mom_probs==max(mom_probs))-1)
    return(mom_probs)
}



