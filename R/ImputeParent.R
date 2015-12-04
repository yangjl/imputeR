#' \code{Impute Parent Genotype}
#' 
#' We have observed parent and observed kids. We want to know the likelihood of the parent's possible genotypes. 
#' This function is to impute the parent's most likely genotype from a progeny array of k kids by 
#' giving a log likelihood threshold.
#' 
#' @param GBS.array A GBS.array object generated from create_array() or sim.array() functions.
#' In the object, a vector (at slot freq) of allele frequencies for each locus will be provided. 
#' It was estimated from the selfed progeny.
#' 
#' @param parents A list of the genotypes of all possible parents, including current parent. 
#' These are observed GBS data, coded with 0, 1 and 2, which is the copy of alternative alleles. 
#' Missing data should be coded with 3.
#' @param imiss Individual missing rate. It used for the error matrices.
#' @param obs_parent An integer referring to which parent is the current parent of interest.
#' @param other_parents A vector of integers referring to which are the second parent of each kid. 
#' Should be same length as number of kids (length(obs_kids))
#' @param obs_kids A list of vectors of Kid's GBS data. Coded with 0, 1 and 2, which is the copy of alternative alleles. 
#' Missing data should be coded with 3.
#' @param hom.error Homozygous error rate, default=0.02.
#' @param het.error Heterozygous error rate, default=0.8.
#' @param oddratio The cutoff used for determining the log likelihood ratio of the highest and the 2nd highest genotypes. 
#' The oddratio = NULL means to report the most likely genotype. Default value sets to 0.6931472, 
#' which is equivalent to most likely genotype being twice as likely as next most likely. 2X as likely
#' @param returnall For function 'parentgeno', returnall=TRUE will return with all the information; returnall=FALSE will return
#' parent's genotype only (either maxlikelihood or, if oddratio is specified, genotypes with low oddratios will be set to missing).
#' @return return a data.frame with all the log likelihoods. Using function \code{parentgeno} to extract parent's genotype. 
#'   
#'   See \url{https://github.com/yangjl/imputeR/blob/master/vignettes/imputeR-vignette.pdf} for more a simulated example and more details
#' @examples
#' test <- sim.array(size.array=50, numloci=10)
#' res <- impute_parent(GBS.array=test)
#' out <- parentgeno(res, oddratio=0.69, returnall=TRUE)
#'
impute_parent <- function(GBS.array, hom.error=0.02, het.error=0.8, imiss=0.5){
        
    ### need to check genotypes
    numloci <- length(GBS.array@gbs_parents[[1]])
    ped <- GBS.array@pedigree
    
    message(sprintf("###>>> Loading a progeny array with [ %s ] GBS loci", numloci))
    message(sprintf("###>>> Of [ %s ] kids, [ %s ] are outcrossed and [ %s ] are selfed", nrow(ped),
                    nrow(subset(ped, p1!=p2)), nrow(subset(ped, p1==p2))))
    ### get probability matrices
    gen_error <- gen_error_mat(hom.error, het.error, imiss)
    probs <- error_mx(hom.error, het.error, imiss)
    
    ### make sfs if not provided?
    p <- GBS.array@freq
    if(is.null(p)){
        message(sprintf("###>>> no allele frequencies provided. Generating random allele frequencies from a neutral SFS"))
        sfs <- getsfs(x = 1:99/100)
        p <- sample(sfs, numloci,replace=TRUE) #get freqs for all loci
    }
    
    #### assign values to impute_one_site function.
    parents <- GBS.array@gbs_parents
    obs_parent <- unique(ped$p1)
    other_parents <- ped$p2
    obs_kids <- GBS.array@gbs_kids
    
    res <- lapply(1:numloci, function(locus){
        impute_one_site(locus, gen_error, p[locus], probs, parents, obs_parent, other_parents, obs_kids)
        })
    geno <- as.data.frame(matrix(unlist(res), ncol=3, byrow=TRUE))
    names(geno) <- c("g0", "g1", "g2")
    
    return(geno)
}

#' @rdname impute_parent
impute_one_site <- function(locus, gen_error, p_locus, probs, parents, obs_parent, other_parents, obs_kids){
    
    obs_parent_probs <- as.numeric()
    for(inferred_parent in 1:3){
        #P(G'obs_parent|G)
        pg_obs <- gen_error[inferred_parent, parents[[obs_parent]][locus]+1] #+1 because obs_parent is 0,1, or 2
    
        #P(G)
        pg <- hw_probs(p_locus)[inferred_parent]
                
        #P(kids|G) sum of logs instead of product
        pkg <- sum(sapply(1:length(obs_kids), function(z){
            ifelse(other_parents[z]==obs_parent,
                log(sum(probs[[inferred_parent]][[inferred_parent]][, obs_kids[[z]][locus]+1])),    
                log(sum(probs[[which.max(sapply(1:3, function(second_parent)  
                    log(hw_probs(p_locus)[second_parent]) + 
                    log(gen_error[second_parent, parents[[other_parents[z]]][locus]+1])+
                    log(sum(probs[[second_parent]][[inferred_parent]][, obs_kids[[z]][locus]+1]))))]][[inferred_parent]][, obs_kids[[z]][locus]+1]))   
            )
        } ))    
        obs_parent_probs[inferred_parent] <- pkg+log(pg_obs)+log(pg)
    }
    return(obs_parent_probs)
}

pkg <- 

#explain this:
#log(sum(probs[[which.max(sapply(1:3, function(second_parent) 
#log(hw_probs(p_locus)[second_parent])+ log(gen_error[second_parent, parents[[other_parents[z]]][locus]+1])+
#log(sum(probs[[second_parent]][[inferred_parent]][, obs_kids[[z]][locus]+1]))))]][[inferred_parent]][, obs_kids[[z]][locus]+1]))   
# 
# we're finding dad that maximizes likelihood of kid z and mom.
# note that this isn't really dad,mom but instead focal parent (which we call mom) and other parents (which we call dads)
# we do this for three dads over variable second_parent
#
# log(hw_probs(p_locus)[second_parent]) <<<< hw prob of a particular dad P(G)
#
# log(gen_error[second_parent, parents[[other_parents[z]]][locus]+1]) 
#<<<<< prob. of observed dad genotype given particular dad (P(G'|G))
#
# log(sum(probs[[second_parent]][[inferred_parent]][, obs_kids[[z]][locus]+1])) 
#<<<<< sum across unknown kid z's genotypes for particular dad (P(G_k|G_mom,G_dad))
# 
# sum the above and you get "blah" (see below)
# which.max(sapply(1:3, function(second_parent) blah )) 
#<<<< gives number of dad (out of 1:3) maximizes likelihood. call this ml_dad
# we then put ml_dad into:
# log(sum(probs[[ml_dad]][[inferred_parent]][, obs_kids[[z]][locus]+1]))   
#<<<< sum across unknown kid's genotype for particular kid z, ml_dad, and inferred_parent (mom) 

#' @rdname impute_parent
#' @param geno A table of genotype likelihoods from impute_parent
#' @param oddratio The cutoff used for determining the log likelihood ratio of the highest and the 2nd highest genotypes. 
#' The oddratio = NULL means to report the most likely genotype. Default value sets to 0.6931472, 
#' which is equivalent to most likely genotype being twice as likely as next most likely. 2X as likely
#' @param returnall If TRUE this will return all the data. 
#' Otherwise it will return the genotypes that pass the specified oddratio; 
#' if no oddratio is specied (oddratio=NULL) then this returns the maximum likelihood genotypes
#' 
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





