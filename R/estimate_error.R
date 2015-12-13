#'
#' \code{Create GBS.array objects} 
#'
#' Create GBS.array objects from Geno4imputeR objects. It will compute population allele frq based on families
#' of selfed kids. Original Geno4imputeR will be updated and returned. A set of GBS.array objects will be output
#' to the directory by chromosomes.
#' 
#' @param Geno4imputeR input object.
#' @param ped a data.frame of pedigree information, must include three columns: proid (progeny id),
#' parent1 (the first parent id), parent2 (the 2nd parent id). 
#' @param outdir path of the output dir.
#' @param maf_cutoff cutoff for MAF, default=0.002.
#' @param lmiss_cutoff cutoff for locus missing rate, default=0.8.
#' @param imiss_cutoff cutoff for individual missing rate, default=0.8.
#' @param size_cutoff Minimum family size required for imputation, default=40.
#' 
#' @return Return updated Geno4imputeR object. 
#' The GBS.array objects will be output to the given directory by family and chr.
#' 
#' Details of the GBS.array objects.
#' Slot1: true_parents, a list of data.frame(hap1, hap2).
#' Slot2: gbs_parents, a list of genotypes. For example, c(1, 2, 2, 0, 3).
#' Slot3: true_kids, a list of data.frame(hap1, hap2).
#' Slot4: gbs_kids, a list of kid genotypes. For example, c(1, 1, 3, 1, 2).
#' Slot5: pedigree, a data.frame (kid, p1, p2). Note, p1 is the focal parent.
#' Slot6: freq, a vector of reference allele freq for all SNPs.
#' 
#'   See \url{https://github.com/yangjl/imputeR/blob/master/vignettes/imputeR-vignette.pdf} for more details.
#'   
#' @examples
#' create_array(Geno4imputeR, ped, outdir="largedata/", 
#' maf_cutoff=0.002, lmiss_cutoff=0.8, imiss_cutoff=0.8, size_cutoff=40)
#' 
estimate_error <- function(geno, ped, self_cutoff, depth_cutoff){
    
    geno <- as.data.frame(geno)
    message(sprintf("###>>> Loaded [ %s ] biallelic loci for [ %s ] plants", nrow(geno), ncol(geno) -3))
    
    ped[,1:3] <- apply(ped[,1:3], 2, as.character)
    if( sum(!unique(c(ped$proid, ped$parent1, ped$parent2)) %in% names(geno)) > 0 ){
        stop("###>>> Some plant ID could not be found in the genotype file !")
    }else{
        pinfo <- pedinfo(ped)
        pinfo <- pinfo[order(pinfo$tot, decreasing=TRUE),]
    }
    
    pinfo <- subset(pinfo, nselfer > self_cutoff)
    #message(sprintf("###>>> Calculating error rate using [ %s ] selfed families ...", nrow(pinfo) )
    
    out <- data.frame()
    for(i in 1:nrow(pinfo)){
        ### get pedigree and idx for the p1 and p2
        pid <- as.character(pinfo$founder[i])
        kid <- subset(ped, parent1 == pid & parent2 == pid)$proid
        subgeno <- geno[, c(pid, kid)]
        
        message(sprintf("###>>> computing [ %s ] selfed families ... ", i))
        err <- check_error(subgeno, depth_cutoff)
        err$fam <- pid
        out <- rbind(out, err)
    }
    
    return(out)
}

#'
check_error <- function(subgeno, depth_cutoff){
    
    totcol <- ncol(subgeno)
    subgeno$count0 <- apply(subgeno[, 2:totcol], 1, function(x) sum(x==0))
    subgeno$count1 <- apply(subgeno[, 2:totcol], 1, function(x) sum(x==1))
    subgeno$count2 <- apply(subgeno[, 2:totcol], 1, function(x) sum(x==2))
    subgeno$count3 <- apply(subgeno[, 2:totcol], 1, function(x) sum(x==3))
    
    out <- data.frame()
    
    sub0 <- subset(subgeno, count0 > depth_cutoff & count1 ==0 & count2 ==0)[,1]
    sub1 <- subset(subgeno, count0 > depth_cutoff & count2 > depth_cutoff)[,1]
    sub2 <- subset(subgeno, count0 + count1 == 0 & count2 > depth_cutoff)[,1]
    
    if(length(sub0) > 0){
        sub0 <- sub0[sub0 != 3]
        err0 <- sum(sub0 != 0)/length(sub0)
    }else{
        err0 <- NA
    }
    
    if(length(sub1) > 0){
        sub1 <- sub1[sub1 != 3]
        err1 <- sum(sub1 != 1)/length(sub1)
    }else{
        err1 <- NA
    }
    if(length(sub2) > 0){
        sub2 <- sub2[sub2 != 3]
        err2 <- sum(sub2 != 2)/length(sub2)
    }else{
        err2 <- NA
    }
    
    return(data.frame(er0 = err0, er1=err1, er2=err2))
    
}

