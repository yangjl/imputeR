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
error_mx <- function(major.error, het.error, minor.error, merr){
    mx <- gen_error_mat(major.error, het.error, minor.error)
    probs <- vector("list",3)
    
    ### remember 4th column is missing data NOT SURE THIS SHOULD BE 1!!!
    
    #AA by AA (1,0,0), Aa (1/2,1/2,0), aa (0,1,0)
    probs[[1]] <- list(cbind(mx*matrix(c(1, 0, 0), nrow = 3, byrow=F, ncol=3), merr),
                       cbind(mx*matrix(c(1/2, 1/2, 0), nrow = 3, byrow=F, ncol=3), merr),
                       cbind(mx*matrix(c(0, 1, 0), nrow = 3, byrow=F, ncol=3), merr) )
    #Aa by AA (1/2,1/2,0), Aa (1/4,1/2,1/4), aa (0,1/2,1/2)
    probs[[2]] <- list(cbind(mx*matrix(c(1/2,1/2,0), nrow = 3, byrow=F, ncol=3),merr),
                       cbind(mx*matrix(c(1/4,1/2,1/4), nrow = 3, byrow=F, ncol=3),merr),
                       cbind( mx*matrix(c(0,1/2,1/2), nrow = 3, byrow=F, ncol=3),merr) )
    #aa by AA (0,1,0), Aa (0,1/2,1/2), aa (1,0,0)
    probs[[3]] <- list(cbind(mx*matrix(c(0,1,0), nrow = 3, byrow=F, ncol=3),merr),
                       cbind(mx*matrix(c(0,1/2,1/2), nrow = 3, byrow=F, ncol=3),merr),
                       cbind(mx*matrix(c(0,0,1), nrow = 3, byrow=F, ncol=3),merr) )
    
    return(probs)
}

#' @rdname error_mx
gen_error_mat <- function(major.error, het.error, minor.error){
    mx <- matrix(c(1-major.error,major.error/2,major.error/2,het.error/2,
                   1-het.error,het.error/2,minor.error/2,minor.error/2,1-minor.error),
                 byrow=T,nrow=3,ncol=3)
    rownames(mx) <- c("g0", "g1", "g2")
    colnames(mx) <- c("ob0", "ob1", "ob2")
    return(mx)
}

#' @rdname error_mx
error_mx2 <- function(major.error, het.error, minor.error){
    mx <- gen_error_mat(major.error, het.error, minor.error)
    probs <- vector("list",3)
    
    ### remember 4th column is missing data NOT SURE THIS SHOULD BE 1!!!
    
    #AA by AA (1,0,0), Aa (1/2,1/2,0), aa (0,1,0)
    probs[[1]] <- list(cbind(mx*matrix(c(1, 0, 0), nrow = 3, byrow=F, ncol=3), 1),
                       cbind(mx*matrix(c(1/2, 1/2, 0), nrow = 3, byrow=F, ncol=3), 1),
                       cbind(mx*matrix(c(0, 1, 0), nrow = 3, byrow=F, ncol=3), 1) )
    #Aa by AA (1/2,1/2,0), Aa (1/4,1/2,1/4), aa (0,1/2,1/2)
    probs[[2]] <- list(cbind(mx*matrix(c(1/2,1/2,0), nrow = 3, byrow=F, ncol=3),1),
                       cbind(mx*matrix(c(1/4,1/2,1/4), nrow = 3, byrow=F, ncol=3),1),
                       cbind( mx*matrix(c(0,1/2,1/2), nrow = 3, byrow=F, ncol=3),1) )
    #aa by AA (0,1,0), Aa (0,1/2,1/2), aa (1,0,0)
    probs[[3]] <- list(cbind(mx*matrix(c(0,1,0), nrow = 3, byrow=F, ncol=3),1),
                       cbind(mx*matrix(c(0,1/2,1/2), nrow = 3, byrow=F, ncol=3),1),
                       cbind(mx*matrix(c(0,0,1), nrow = 3, byrow=F, ncol=3),1) )
    
    return(probs)
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
hw_probs <- function(x){ return(c((1-x)^2,2*x*(1-x),x^2))}
############################################################################

kids_errs <- function(simk, imputek){
    #
    out <- data.frame()
    for(i in 1:length(simk)){
        
        imputeki <- imputek[[i]][[3]]
        simki <- simk[[i]][[2]][imputeki$idx, ]
        
        names(simki) <- c("shap1", "shap2", "obs")
        comb0 <- cbind(simki, imputeki)
        comb <- subset(comb0, k1 != 3)
        err = 0
        for(j in unique(comb$chunk)){
            sub <- subset(comb, chunk == j)
            idx1 <- which.max(c(cor(sub$k1, sub$shap1), cor(sub$k1, sub$shap2)) )
            err1 <- sum(sub$k1 != sub[, idx1])
            idx2 <- which.max(c(cor(sub$k2, sub$shap1), cor(sub$k2, sub$shap2)) )
            err2 <- sum(sub$k2 != sub[, idx2])
            err <- err + err1 + err2
        }
        tem <- data.frame(kid=i, chunks=length(unique(comb$chunk)), err=err, tot=nrow(simki), miss= 1-nrow(comb)/nrow(comb0))
        tem$geno <- sum(comb$shap1+comb$shap2 != comb$k1+comb$k2)
        out <- rbind(out, tem)
    }
    out$phaserate <- out$err/out$tot
    out$genorate <- out$geno/out$tot
    return(out)
}

#### phasing results
mom_phasing_error <- function(phasemom, sim){
    names(phasemom)[3:4] <- c("ihap1", "ihap2")
    
    #### error estimation
    outall <- cbind(phasemom, sim[[1]][phasemom$idx, ])
    err = 0
    for(i in unique(outall$chunk)){
        chunk <- subset(outall, chunk == i)
        idx1 <- which.max(c(cor(chunk$ihap1, chunk$hap1), cor(chunk$ihap1, chunk$hap2)) )
        err1 <- sum(chunk$ihap1 != chunk[, 4+idx1])
        idx2 <- which.max(c(cor(chunk$ihap2, chunk$hap1), cor(chunk$ihap2, chunk$hap2)) )
        err2 <- sum(chunk$ihap2 != chunk[, 4+idx2])
        err <- err + err1 + err2
    }
    out <- data.frame(chunk=length(unique(outall$chunk)), err=err, loci=2*nrow(outall))
    out$rate <- out$err/out$loci
    return(out)
}
#### dad phasing results
dad_phasing_error <- function(newdad, simdad){
    names(newdad)[3:4] <- c("ihap1", "ihap2")
    
    #### error estimation
    outall <- cbind(newdad, simdad[newdad$idx, ])
    err = 0
    for(i in unique(outall$chunk)){
        chunk <- subset(outall, chunk == i)
        idx1 <- which.max(c(cor(chunk$ihap1, chunk$hap1), cor(chunk$ihap1, chunk$hap2)) )
        err1 <- sum(chunk$ihap1 != chunk[, 4+idx1])
        idx2 <- which.max(c(cor(chunk$ihap2, chunk$hap1), cor(chunk$ihap2, chunk$hap2)) )
        err2 <- sum(chunk$ihap2 != chunk[, 4+idx2])
        err <- err + err1 + err2
    }
    out <- data.frame(chunk=length(unique(outall$chunk)), err=err, loci=2*nrow(outall))
    out$rate <- out$err/out$loci
    return(out)
}

#' \code{Method for GBS.array} 
#'
#' Simulate imputed parents. All the gbs_parents are true genotypes.
#' 
#' @param GBS.array Input a GBS.array object.
#' @param phased.parents Whether use the phased non-focal parents.
#' 
get_true_GBS <- function(GBS.array, phased.parents=TRUE){
    
    ped <- GBS.array@pedigree
    if(length(unique(ped$p1)) != 1){
        stop("### more than one focal parent!!!")
    }
    
    if(phased.parents){
        idx <- unique(GBS.array@pedigree$p1)
        pidx <- 1:length(GBS.array@true_parents)
        pidx <- pidx[pidx != idx]
        
        for(i in pidx){
            tem <- GBS.array@true_parents[[i]]
            tem$idx <- 1:nrow(tem)
            tem$chunk <- 1
            GBS.array@gbs_parents[[i]] <- tem
        }
        
        true_p <- GBS.array@true_parents[[idx]]
        GBS.array@gbs_parents[[idx]] <- true_p$hap1 + true_p$hap2
        
    }else{
        for(pidx in 1:length(GBS.array@true_parents)){
            true_p <- GBS.array@true_parents[[pidx]]
            GBS.array@gbs_parents[[pidx]] <- true_p$hap1 + true_p$hap2
        }
    }
    return(GBS.array)
}

#' \code{Method for GBS.array} 
#'
#' Simulate phased parents. All the gbs_parents are perfectly phased.
#' 
#' @param GBS.array Input a GBS.array object.
#' @param maxchunk The max number of chunks.
#' 
get_phased <- function(GBS.array, maxchunk){
    
    ped <- GBS.array@pedigree
    if(length(unique(ped$p1)) != 1){
        stop("### more than one focal parent!!!")
    }
    
    for(pidx in 1:length(GBS.array@true_parents)){
        true_p <- GBS.array@true_parents[[pidx]]
        c <- sample(1:maxchunk, 1)
        idx <- sample(1:c, nrow(true_p), replace=TRUE)
        true_p$chunk <- sort(idx)
        true_p$idx <- 1:nrow(true_p)
        GBS.array@gbs_parents[[pidx]] <- true_p 
    }
    return(GBS.array)
}