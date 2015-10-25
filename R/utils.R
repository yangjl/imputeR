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

    #row/col names <- currently not used, need to fix for clarity
    rown <- c("g0", "g1", "g2")
    coln <- c("ob0", "ob1", "ob2", "obN")
    
    ### remember 4th column is missing data NOT SURE THIS SHOULD BE 1!!!
    
    #AA by AA (1,0,0), Aa (1/2,1/2,0), aa (0,1,0)
    probs[[1]] <- list(cbind(mx*matrix(c(1, 0, 0), nrow = 3,byrow=F,ncol=3),1/3),cbind(mx*matrix(c(1/2, 1/2, 0), nrow = 3,byrow=F,ncol=3),1/3),cbind(mx*matrix(c(0, 1, 0), nrow = 3,byrow=F,ncol=3),1/3)   )
    #Aa by AA (1/2,1/2,0), Aa (1/4,1/2,1/4), aa (0,1/2,1/2)
    probs[[2]] <- list(cbind(mx*matrix(c(1/2,1/2,0), nrow = 3,byrow=F,ncol=3),1/3),cbind(mx*matrix(c(1/4,1/2,1/4), nrow = 3,byrow=F,ncol=3),1/3),cbind(mx*matrix(c(0,1/2,1/2), nrow = 3,byrow=F,ncol=3),1/3)   )
    #aa by AA (0,1,0), Aa (0,1/2,1/2), aa (1,0,0)
    probs[[3]] <- list(cbind(mx*matrix(c(0,1,0), nrow = 3,byrow=F,ncol=3),1/3),cbind(mx*matrix(c(0,1/2,1/2), nrow = 3,byrow=F,ncol=3),1/3),cbind(mx*matrix(c(0,0,1), nrow = 3,byrow=F,ncol=3),1/3)   )
    
    return(probs)
}

#' @rdname Utils
gen_error_mat <- function(hom.error, het.error){
    mx <- matrix(c(1-hom.error,hom.error/2,hom.error/2,het.error/2,
                   1-het.error,het.error/2,hom.error/2,hom.error/2,1-hom.error),
                 byrow=T,nrow=3,ncol=3)
    mx <- cbind(mx, 1)
    rownames(mx) <- c("g0", "g1", "g2")
    colnames(mx) <- c("ob0", "ob1", "ob2", "obN")
    return(mx)
}

#' @rdname Utils 
# Create random haplotype with sfs
ran.hap <- function(numloci,p){
    sapply(1:numloci,function(x) rbinom(1,1,p[x]))
}

#' @rdname Utils 
# Add error to diploid
add_error<-function(diploid,hom.error,het.error){
    hets_with_error=sample(which(diploid==1),round(het.error*length(which(diploid==1))))
    hom0_with_error=sample(which(diploid==0),round(hom.error*length(which(diploid==0))))
    hom2_with_error=sample(which(diploid==2),round(hom.error*length(which(diploid==2))))
    diploid=replace(diploid,hets_with_error,sample(c(0,2),length(hets_with_error),replace=T)  )
    ### error rate from (hom => het) == (hom1 => hom2)
    diploid=replace(diploid,hom0_with_error,sample(c(1,2),length(hom0_with_error),replace=T)  )
    diploid=replace(diploid,hom2_with_error,sample(c(1,0),length(hom2_with_error),replace=T)  )
    return(diploid)
}

#' @rdname Utils 
#Copy mom to kids with recombination
copy.mom <- function(mom, co_mean){ 
    co=rpois(1,co_mean) #crossovers
    numloci=length(mom[[1]])
    recp=c(1,sort(round(runif(co, min=2, max=numloci-1))), numloci+1) #position   
    chrom=rbinom(1,1,.5)+1
    kpiece=as.numeric()
    hap <- c()
    for(r in 1:(length(recp)-1)){
        kpiece=c(kpiece,mom[[chrom]][recp[r]:(recp[r+1]-1)]) #copy 1->rec from mom
        hap <- c(hap, chrom)
        chrom=ifelse(chrom==1,2,1)  
    }     
    return(list(kpiece, data.frame(hap=hap, start=recp[-length(recp)], end=recp[-1]) ))
}

#' @rdname Utils 
# add missing
missing.idx <- function(nloci, imiss){
    #hist(rbeta(10000, 2, 2))
    if(imiss >= 1){
        m <- rbeta(1, 2, 2)
    }else{
        m <- imiss
    }
    
    ml <- sort(sample(1:nloci, size=round(m*nloci)))
    return(ml)
}


#' @rdname Utils 
# Returns a list of true [[1]] and observed [[2]] kid
kid <- function(mom, dad, het.error, hom.error, rec=1.5, imiss=0.3, misscode=3){
    if(rec==0){
        k1=mom[[rbinom(1,1,.5)+1]]
        k2=dad[[rbinom(1,1,.5)+1]]
    } else{
        k1=copy.mom(mom,rec) # list
        k2=copy.mom(dad,rec)
    }
    true_kid=k1[[1]] + k2[[1]]
    #return(list(true_kid,obs_kid))
    
    obs_kid <- add_error(true_kid, hom.error, het.error)
    if(imiss > 0){
        idx <- missing.idx(length(true_kid), imiss)
        obs_kid <- replace(obs_kid, idx, misscode)
    }
    
    info <- list(k1[[2]], k2[[2]])
    simk <- data.frame(hap1=k1[[1]], hap2=k2[[1]], obs= obs_kid )
    return(list(info, simk))
}
