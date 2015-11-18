#' Return error matrices
#'
#' @param hom.error homozygous error rate, default=0.02.
#' @param het.error heterozygous error rate, default=0.8.
#' @return \code{error_mx} returns a list of three matrices.
#' In this matrix,
#' row1 is true_gen 00, row2 is true_gen 01, row3 is true_gen 11.
#' cols 1-3 are obs. genotype (00,01,11) and last col4 is the missing data.
#' 1. 00x01
#' 2. 01x01
#' 3. 11x01
#' @examples
#' #Genotype error matrix
#' gen_error_mat(hom.error=0.02, het.error=0.8)
#' 
#' #Genotype error by consideringMendelian segregation rate
#' error_mx(hom.error=0.02, het.error=0.8)
#'
error_mx <- function(hom.error, het.error, imiss){
    mx <- gen_error_mat(hom.error, het.error, imiss)
    probs <- vector("list",3)
    
    #row/col names <- currently not used, need to fix for clarity
    rown <- c("g0", "g1", "g2")
    coln <- c("ob0", "ob1", "ob2", "obN")
    
    ### remember 4th column is missing data NOT SURE THIS SHOULD BE 1!!!
    
    #AA by AA (1,0,0), Aa (1/2,1/2,0), aa (0,1,0)
    probs[[1]] <- list(mx*matrix(c(1, 0, 0), nrow = 3, byrow=F, ncol=4),
                       mx*matrix(c(1/2, 1/2, 0), nrow = 3, byrow=F, ncol=4),
                       mx*matrix(c(0, 1, 0), nrow = 3, byrow=F, ncol=4) )
    #Aa by AA (1/2,1/2,0), Aa (1/4,1/2,1/4), aa (0,1/2,1/2)
    probs[[2]] <- list(mx*matrix(c(1/2,1/2,0), nrow = 3, byrow=F, ncol=4),
                       mx*matrix(c(1/4,1/2,1/4), nrow = 3, byrow=F, ncol=4),
                       mx*matrix(c(0,1/2,1/2), nrow = 3, byrow=F, ncol=4) )
    #aa by AA (0,1,0), Aa (0,1/2,1/2), aa (1,0,0)
    probs[[3]] <- list(mx*matrix(c(0,1,0), nrow = 3, byrow=F, ncol=4),
                       mx*matrix(c(0,1/2,1/2), nrow = 3, byrow=F, ncol=4),
                       mx*matrix(c(0,0,1), nrow = 3, byrow=F, ncol=4) )
    
    return(probs)
}

#' @rdname error_mx
gen_error_mat <- function(hom.error, het.error, imiss){
    mx <- matrix(c(1-hom.error,hom.error/2,hom.error/2,het.error/2,
                   1-het.error,het.error/2,hom.error/2,hom.error/2,1-hom.error),
                 byrow=T,nrow=3,ncol=3)
    mx <- cbind(mx*(1-imiss), imiss)
    rownames(mx) <- c("g0", "g1", "g2")
    colnames(mx) <- c("ob0", "ob1", "ob2", "obN")
    return(mx)
}


#' setup the neutral SFS
#'
#' @param x A vector of freq bins. Default 0.01-0.99.
#' @return Return a vector of neutral SFS (N=51734)
#' @examples
#' sfs <- getsfs(x=1:99/100)
#' hist(sfs)
#'
getsfs <- function(x=1:99/100){
    freq <- 1/x
    sfs <- as.numeric(); 
    for(i in 1:99){
        sfs <- c(sfs,rep(x[i],100*freq[i]))
    }
    return(sfs)
}

#' Return HW probs
#' 
#' @param x Allele freq p.
#' @return Return a vector of p^2, 2pq and q^2,
#' @examples
#' hw_probs(0.1)
#'
hw_probs <- function(x){ return(c(x^2,2*x*(1-x),(1-x)^2))}
############################################################################

