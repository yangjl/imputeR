#' setup the neutral SFS
#'
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

# Return HW probs
hw_probs <- function(x){ return(c(x^2,2*x*(1-x),(1-x)^2))}
############################################################################

#' Return error matrices
#'
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
error_mx <- function(hom.error, het.error){
    
    mx <- gen_error_mat(hom.error=hom.error, het.error=het.error)[, -4]
    probs <- vector("list",3)

    #row/col names
    rown <- c("g0", "g1", "g2")
    coln <- c("ob0", "ob1", "ob2", "obN")
    
    ### 1: 00x01; 2: 01x01; 3: 11x01
    probs[[1]] <-  cbind(mx*matrix(c(1/2, 1/2, 0), nrow = 3,byrow=F,ncol=3),1); row.names(probs[[1]])=rown; colnames(probs[[1]])=coln
    probs[[2]] <- cbind(mx*matrix(c(1/4, 1/2, 1/4), nrow = 3,byrow=F,ncol=3),1); row.names(probs[[2]])=rown; colnames(probs[[2]])=coln
    probs[[3]] <- cbind(mx*matrix(c(0, 1/2, 1/2), nrow = 3,byrow=F,ncol=3),1); row.names(probs[[3]])=rown; colnames(probs[[3]])=coln
    return(probs)
}

#' @rdname error_mx
gen_error_mat <- function(hom.error, het.error){
    mx <- matrix(c(1-hom.error,hom.error/2,hom.error/2,het.error/2,
                   1-het.error,het.error/2,hom.error/2,hom.error/2,1-hom.error),
                 byrow=T,nrow=3,ncol=3)
    mx <- cbind(mx, 1)
    rownames(mx) <- c("g0", "g1", "g2")
    colnames(mx) <- c("ob0", "ob1", "ob2", "obN")
    return(mx)
}

