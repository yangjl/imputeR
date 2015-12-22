#'
#' \code{Estimate genotyping error} 
#'
#' Estimate error from data and pedigree.
#' 
#' @param geno input object.
#' @param ped a data.frame of pedigree information, must include three columns: proid (progeny id),
#' parent1 (the first parent id), parent2 (the 2nd parent id). 
#' @param imiss_cutoff cutoff for individual missing rate, default=0.8.
#' @param size_cutoff Minimum family size required for imputation, default=40.
#' 
#' @return Return error matrix.
#'   See \url{https://github.com/yangjl/imputeR/blob/master/vignettes/imputeR-vignette.pdf} for more details.
#'   
#' @examples
#' library(data.table, lib="~/bin/Rlib/")
#' ped <- read.table("data/parentage_info.txt", header =TRUE)
#' ped[, 1:3] <- apply(ped[, 1:3], 2, as.character)
#' geno <- fread("largedata/lcache/teo_recoded.txt")
#' estimate_error(geno, ped, self_cutoff=30, depth_cutoff=10)
#' 
estimate_error <- function(geno, ped, self_cutoff, depth_cutoff, est_kids=FALSE){
    
    if( sum(!unique(c(ped$proid, ped$parent1, ped$parent2)) %in% names(geno)) > 0 ){
        stop("###>>> Some plant ID could not be found in the genotype file !")
    }else{
        #geno <- as.data.frame(geno)
        message(sprintf("###>>> Loaded [ %s ] biallelic loci for [ %s ] plants", nrow(geno), nrow(ped)))

        pinfo <- pedinfo(ped)
        pinfo <- pinfo[order(pinfo$tot, decreasing=TRUE),]
    }
    
    pinfo <- subset(pinfo, nselfer > self_cutoff)
    message(sprintf("###>>> Calculating error rate using [ %s ] selfed families ...", nrow(pinfo) ))
    
    out1 <- data.frame()
    for(i in 1:nrow(pinfo)){
        ### get pedigree and idx for the p1 and p2
        pid <- as.character(pinfo$founder[i])
        kid <- subset(ped, parent1 == pid & parent2 == pid)$proid
        subgeno <- geno[, c(pid, kid)]
        
        message(sprintf("###>>> computing [ %s ] selfed families ... ", i))
        err <- check_error(subgeno, depth_cutoff)
        err$fam <- pid
        out1 <- rbind(out1, err)
    }
    
    if(est_kids){
        out2 <- kid_het_err(geno, depth_cutoff, ped, pinfo, verbose=TRUE)
        return(list(out1, out2))
    }else{
        return(out1)
    }
    
}

#' @rdname estimate_error
kid_het_err <- function(geno, depth_cutoff, ped, pinfo, verbose){
    
    if(verbose) message(sprintf("###>>> start to compute het err."))
    #### matrix of genotype counts using selfed kids
    count_mx <- list()
    for(i in 1:nrow(pinfo)){
        ### get pedigree and idx for the p1 and p2
        pid <- as.character(pinfo$founder[i])
        kid <- subset(ped, parent1 == pid & parent2 == pid)$proid
        subgeno <- geno[, c(pid, kid)]
        count_mx[[pid]] <- count_geno(subgeno)[, (ncol(subgeno)+1):(ncol(subgeno)+3)]
    }
    
    ### others
    ocped <- subset(ped, parent1 != parent2 & (parent1 %in% pinfo$founder & parent2 %in% pinfo$founder))
    
    out <- data.frame()
    for(j in 1:nrow(ocped)){
        if(verbose) message(sprintf("###>>> processing the [ %s ] kids", j))
        p1_idx0 <- subset(count_mx[[ocped$parent1[j]]], count0 > depth_cutoff & count1 + count2 ==0)
        p2_idx0 <- subset(count_mx[[ocped$parent2[j]]], count0 > depth_cutoff & count1 + count2 ==0)
        
        p1_idx2 <- subset(count_mx[[ocped$parent1[j]]], count0 + count1 ==0 & count2 > depth_cutoff)
        p2_idx2 <- subset(count_mx[[ocped$parent2[j]]], count0 + count1 ==0 & count2 > depth_cutoff)
        
        ## het sites
        h1 <- merge(p1_idx0, p2_idx2, by="row.names")
        h2 <- merge(p1_idx2, p2_idx0, by="row.names")
        hetgeno <- geno[c(h1$Row.names, h2$Row.names) , ocped$proid[j]]
        hetgeno <- hetgeno[hetgeno != 3]
        
        #### major
        maj <- merge(p1_idx0, p2_idx0, by="row.names")
        majgeno <- geno[maj$Row.names, ocped$proid[j]]
        majgeno <- majgeno[majgeno != 3]
        
        #### minor
        mnr <- merge(p1_idx2, p2_idx2, by="row.names")
        mnrgeno <- geno[mnr$Row.names, ocped$proid[j]]
        mnrgeno <- mnrgeno[mnrgeno != 3]
        
        tem <- data.frame(kid=ocped$proid[j], 
                          nhet=length(hetgeno), kerr10=sum(hetgeno==0)/length(hetgeno), kerr12=sum(hetgeno==2)/length(hetgeno),
                          nmaj=length(majgeno), kerr01=sum(majgeno==1)/length(majgeno), kerr02=sum(hetgeno==2)/length(majgeno),
                          nmnr=length(mnrgeno), kerr20=sum(mnrgeno==0)/length(mnrgeno), kerr21=sum(mnrgeno==1)/length(mnrgeno))
        out <- rbind(out, tem)
    }
    
    return(out)
}

#' @rdname estimate_error
count_geno <- function(subgeno){
    
    totcol <- ncol(subgeno)
    subgeno$count0 <- apply(subgeno[, 2:totcol], 1, function(x) sum(x==0))
    subgeno$count1 <- apply(subgeno[, 2:totcol], 1, function(x) sum(x==1))
    subgeno$count2 <- apply(subgeno[, 2:totcol], 1, function(x) sum(x==2))
    #subgeno$count3 <- apply(subgeno[, 2:totcol], 1, function(x) sum(x==3))
    return(subgeno)
}

#' @rdname estimate_error
check_error <- function(subgeno, depth_cutoff){
    
    totcol <- ncol(subgeno)
    subgeno$count0 <- apply(subgeno[, 2:totcol], 1, function(x) sum(x==0))
    subgeno$count1 <- apply(subgeno[, 2:totcol], 1, function(x) sum(x==1))
    subgeno$count2 <- apply(subgeno[, 2:totcol], 1, function(x) sum(x==2))
    subgeno$count3 <- apply(subgeno[, 2:totcol], 1, function(x) sum(x==3))
    
    sub0 <- subset(subgeno, count0 > depth_cutoff & count1 ==0 & count2 ==0)[,1]
    sub1 <- subset(subgeno, count0 > depth_cutoff & count2 > depth_cutoff)[,1]
    sub2 <- subset(subgeno, count0 + count1 == 0 & count2 > depth_cutoff)[,1]
    
    #### founder error
    if(length(sub0) > 0){
        sub0 <- sub0[sub0 != 3]
        err0 <- sum(sub0 != 0)/length(sub0)
        err01 <- sum(sub0 == 1)/length(sub0)
        err02 <- sum(sub0 == 2)/length(sub0)
    }else{
        err0 <- err01 <- err02 <- NA
    }
    
    if(length(sub1) > 0){
        sub1 <- sub1[sub1 != 3]
        err1 <- sum(sub1 != 1)/length(sub1)
        err10 <- sum(sub1 == 0)/length(sub1)
        err12 <- sum(sub1 == 2)/length(sub1)
    }else{
        err1 <- err10 <- err12 <- NA
    }
    if(length(sub2) > 0){
        sub2 <- sub2[sub2 != 3]
        err2 <- sum(sub2 != 2)/length(sub2)
        err20 <- sum(sub2 == 0)/length(sub2)
        err21 <- sum(sub2 == 1)/length(sub2)
    }else{
        err2 <- err20 <- err21 <- NA
    }
    
    #### kid err:
    kid0 <- subset(subgeno, subgeno[,1] == 0)
    k01 <- sum(kid0$count1)/((totcol-1)*nrow(kid0) - sum(kid0$count3))
    k02 <- sum(kid0$count2)/((totcol-1)*nrow(kid0) - sum(kid0$count3))
    
    kid2 <- subset(subgeno, subgeno[,1] == 2)
    k21 <- sum(kid2$count1)/((totcol-1)*nrow(kid2) - sum(kid2$count3))
    k20 <- sum(kid2$count0)/((totcol-1)*nrow(kid2) - sum(kid2$count3))
    
    return(data.frame(er0 = err0, er01=err01, er02=err02,
                      er1=err1, er10=err10, er12=err12,
                      er2=err2, er20=err20, er21=err21,
                      k01=k01, k02=k02, k21=k21, k20=k20))
    
}

