#'
#' \code{Create imputeR object. } 
#'
#' S4 method for TasselHDF5 object. Recode it into `0, 1, 2, 3` format. 
#' In the genotype file, `0, 1, 2` indicate the copy of alternative allele,
#' `3` indicate missing data. This function also calculate individual and locus missing rate.
#' 
#' @param inpob A S4 method for TasselHDF5 object.
#' @return return "Geno4imputeR" object with three slots: "genomx", "info", "imiss". 
#' "genomx", genotype matrix. 
#' "info", loci information, including SNP ID (snpid), chr, position (pos), 
#' loci missing rate (lmiss), minor allele freq (maf). 
#' "imiss", individual plant missing rate.
#'  
#'   See \url{https://github.com/yangjl/imputeR/blob/master/vignettes/imputeR-vignette.pdf} for more details.
#'   
#' @examples
#' # note this function requires at least 100G memeory to load teo.h5 file.
#' obj <- loading_hdf5(h5file="largedata/teo.h5")
#' 
imputeRob <-  function(teo){
    pos <- as.data.frame(granges(teo))
    alt <- sapply(teo@alt, function(x) TASSELL_ALLELES[x+1L])
    info <- data.frame(snpid=rownames(geno(teo)), ref=ref(teo), alt=alt)
    info <- merge(info, pos, by.x="snpid", by.y="row.names")
    
    ### get genotype matrix
    genos <- geno(teo)
    
    ### calculate missing rates
    message("calculating missing rates ...", appendLF=FALSE)
    imiss <- apply(genos, 2, function(x) return(sum(is.na(x))/nrow(genos)))
    imiss <- data.frame(imiss)
    
    lmiss <- apply(genos, 1, function(x) return(sum(is.na(x))/ncol(genos)))
    lmiss <- data.frame(lmiss)
    message("done.")
    
    ### calculate maf
    message("calculating minor allele frq (MAF) ...", appendLF=FALSE)
    maf <- apply(genos, 1, function(x){
        x <- x[!is.na(x)]
        c0 <- sum(x == 0)
        c1 <- sum(x == 1)
        c2 <- sum(x == 2)
        return(min(c(2*c0 +c1, c1 + 2*c2))/(2*(c0 + c1 + c2)) )
    })
    message("done.")
    
    ### turn missing data "NA" into 3 and return the results
    message(sprintf("preparing Geno4imputeR object for [%s] plants and [%s] SNPs ...", 
                    ncol(genos), nrow(genos)))
    maf <- data.frame(maf)
    genos[is.na(genos)] <- 3
    info <- merge(info, lmiss, by.x="snpid", by.y="row.names")
    info <- merge(info, maf, by.x="snpid", by.y="row.names")
    info <- info[, c("snpid", "seqnames", "start", "ref", "alt", "lmiss", "maf")]
    names(info)[1:3] <- c("snpid", "chr", "pos")
    
    obj <- new("Geno4imputeR",
               genomx = genos,
               info = info,
               imiss=imiss)
    message("ALL DONE.")
    return(obj)
}
