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