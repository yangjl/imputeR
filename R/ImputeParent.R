#'
#' \code{Impute parent's genotype using parent and kids' GBS data. } 
#'
#' We have observed parent and observed kids. We want to know the likelihood of the parent's possible genotypes. 
#' This function is to impute the parent's most likely genotype from a progeny array of k kids by 
#' giving a log likelihood threshold.
#'
#' @param parents A list of the genotypes of all possible parents, including current parent. These are observed GBS data, coded with 0, 1 and 2, which is the copy of alternative alleles. 
#' Missing data should be coded with 3.
#' @param obs_parent An integer referring to which parent is the current parent of interest.
#' @param other_parents A vector of integers referring to which are the second parent of each kid. 
#' Should be same length as number of kids (length(obs_kids))
#' @param obs_kids A list of vectors of Kid's GBS data. Coded with 0, 1 and 2, which is the copy of alternative alleles. 
#' Missing data should be coded with 3.
#' @param hom.error Homozygous error rate, default=0.02.
#' @param het.error Heterozygous error rate, default=0.8.
#' @param p vector of allele frequencies at each locus. should be estimated from set of all parents.
#' WARNING: if not supplied, a random neutral SFS is used to generate p at each locus.
#' @param oddratio The cutoff used for determining the log likelihood ratio of the highest and the 2nd highest genotypes. 
#' The oddratio = NULL means to report the most likely genotype. Default value sets to 0.6931472, 
#' which is equivalent to most likely genotype being twice as likely as next most likely. 2X as likely
#' @param returnall For function 'parentgeno', returnall=TRUE will return with all the information; returnall=FALSE will return
#' parent's genotype only (either maxlikelihood or, if oddratio is specified, genotypes with low oddratios will be set to missing).
#' @return return a data.frame with all the log likelihoods. Using function \code{parentgeno} to extract parent's genotype. 
#'   
#'   See \url{https://github.com/yangjl/imputeR/blob/master/vignettes/imputeR-vignette.pdf} for more details.
#' @examples
#' obs_parent <- c(0, 0, 0)
#' obs_kids <- list(c(0, 0, 0), c(0, 0, 0), c(0, 0, 0), c(0, 0, 0), c(0, 0, 0), c(0, 0, 0),
#' c(0, 0, 0), c(0, 0, 0), c(0, 0, 0), c(0, 0, 0), c(0, 0, 0), c(0, 0, 0))
#'
#' impute_parent(obs_parent, obs_kids, hom.error=0.02, het.error=0.8)
#' 
#' parentgeno(geno, oddratio=0.5, returnall=FALSE)
#' 
#' 
impute_parent <- function( parents, obs_parent, other_parents, obs_kids, hom.error=0.02, het.error=0.8,p=NULL){
    
    ### need to check genotypes
    numloci <- length(parents[[obs_parent]])
    
    message(sprintf("###>>> impute parent's genotype using [ %s ] kids ...", length(obs_kids)))
    ### get probability matrices
    gen_error <- gen_error_mat(hom.error, het.error)
    probs <- error_mx(hom.error, het.error)
    
    ### make sfs if not provided?
    if(is.null(p)){
        message(sprintf("###>>> no allele frequencies provided. generating random allele frequencies from a neutral SFS"))
        sfs <- getsfs(x = 1:99/100)
        p <- sample(sfs, numloci,replace=TRUE) #get freqs for all loci
    }
    
    res <- lapply(1:length(parents[[obs_parent]]), function(locus) impute_one_site(locus, gen_error, p[locus], probs, parents, obs_parent, other_parents, obs_kids))
    geno <- as.data.frame(matrix(unlist(res), ncol=3, byrow=TRUE))
    names(geno) <- c("g0", "g1", "g2")
    
    return(geno)
}

#' @rdname impute_parent
impute_one_site <- function(locus, gen_error, p_locus, probs, parents, obs_parent, other_parents, obs_kids)){
    obs_parent_probs <- as.numeric()
    for(inferred_parent in 1:3){
        #P(G'obs_parent|G)
        pg_obs <- gen_error[inferred_parent, parents[[obs_parent]][locus]+1] #+1 because obs_parent is 0,1, or 2
    
        #P(G)
        pg <- hw_probs(p_locus)[inferred_parent]
                
        #P(kids|G) sum of logs instead of product
        pkg <- sum(sapply(1:length(obs_kids), function(z){                
                
            #find ML dad
            s_parent_probs<-as.numeric()
            for(second_parent in 1:3){
                ps <- hw_probs(p_locus)[second_parent] #hw probs of given second parent genotype
                ps_obs <- gen_error[second_parent, parents[[other_parents[z]]][locus]+1] #+1 because obs_parent is 0,1, or 2
                s_k=sum(probs[[second_parent]][[inferred_parent]][, obs_kids[[z]][locus]+1])
                s_parent_probs[second_parent]=log(s_k)+log(ps)+log(ps_obs)
            }
            #if two dads equal, randomly picks
            ml_dad=ifelse(length(which(s_parent_probs==max(s_parent_probs)))>1,sample(which(s_parent_probs==max(s_parent_probs)),1),which(s_parent_probs==max(s_parent_probs)))                   
            log(sum(probs[[ml_dad]][[inferred_parent]][, obs_kids[[z]][locus]+1]))   
        } ))    
        obs_parent_probs[inferred_parent] <- pkg+log(pg_obs)+log(pg)
    }
    return(obs_parent_probs)
}


#' @rdname impute_parent
parentgeno <- function(geno, oddratio=0.6931472, returnall=TRUE){ 
    geno$OR <- apply(geno, 1, function(v){
        n <- length(v)
        return(max(v) - sort(v, partial=n-1)[n-1])  
    })
    geno$gmax <- apply(geno[, 1:3], 1, function(v){
        return(which.max(v)-1)  
    })
    geno$gor <- 3 # 3 is missing data
    geno[geno$OR > oddratio, ]$gor <- geno[geno$OR > oddratio, ]$gmax
    
    if(returnall){
        return(geno)
    }else{
        if(is.null(oddratio)){
            return(geno$gmax)   
        }else{
            return(geno$gor)
        } 
    }
}





