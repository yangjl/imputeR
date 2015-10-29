#'
#' \code{Load HDF5 file and recode the genotype matrix. } 
#'
#' Load HDF5 file and recode into `0, 1, 2, 3` format. In the genotype file, `0, 1, 2` indicate the copy of alternative allele,
#' `3` indicate missing data. This function also calculate individual and locus missing rate.
#' 
#' @param h5file h5 file from tassel. For example, "teo.h5", 
#' @param save.file RData object to save. For example, "out.RData".
#' @return return "Geno4imputeR" object with three slots: "genomx", "info", "imiss". 
#' "genomx", genotype matrix. "info", loci information, including missing
#' rate, maf, snp position, ref and alt alleles, etc. "imiss", individual plant missing rate.
#'  
#'   See \url{https://github.com/yangjl/imputeR/blob/master/vignettes/imputeR-vignette.pdf} for more details.
#'   
#' @examples
#' # note this function requires at least 100G memeory to load teo.h5 file.
#' obj <- loading_hdf5(h5file="largedata/teo.h5")
#' 
loading_hdf5 <- function(h5file="largedata/teo.h5"){
   
    #### load h5file
    teo <- initTasselHDF5(h5file, version="5")
    teo <- loadBiallelicGenotypes(teo, verbose = TRUE)
    
    pos <- as.data.frame(granges(teo))
    alt <- sapply(teo@alt, function(x) TASSELL_ALLELES[x+1L])
    info <- data.frame(snpid=rownames(geno(teo)), ref=ref(teo), alt=alt)
    info <- merge(info, pos, by.x="snpid", by.y="row.names")
    
    ### get genotype matrix
    genos <- geno(teo)
    
    ### calculate missing rates
    imiss <- apply(genos, 2, function(x) return(sum(is.na(x))/598043))
    imiss <- data.frame(imiss)
    
    lmiss <- apply(genos, 1, function(x) return(sum(is.na(x))/4875))
    lmiss <- data.frame(lmiss)
    
    ### calculate maf
    maf <- apply(genos, 1, function(x){
        x <- x[!is.na(x)]
        c0 <- sum(x == 0)
        c1 <- sum(x == 1)
        c2 <- sum(x == 2)
        return(min(c(2*c0 +c1, c1 + 2*c2))/(2*(c0 + c1 + c2)) )
    })
    
    maf <- data.frame(maf)
    genos[is.na(genos)] <- 3
    
    ## return the results
    info <- merge(info, lmiss, by.x="snpid", by.y="row.names")
    info <- merge(info, maf, by.x="snpid", by.y="row.names")
    
    obj <- new("Geno4imputeR", filename=file,
               genomx = genos,
               info = info,
               imiss=imiss)
    return(obj)
}
