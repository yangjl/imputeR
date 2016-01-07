#' \code{Phasing and Imputing Kid Genotypes}
#'
#' Calculate the maximum likelihood of the inherited haplotypes. And then find the minimum path of recombinations.
#' At the potential recombination sites, find the maximum likelihood of the breaking points.
#'
#' @param GBS.array A GBS.array object with phased parents. 
#' It can be generated from \code{phase_parent} or\code{sim.array} functions.
#' @param winsize The size of window for determining kid haplotype. Default=100. 
#' @param verbose Writing verbose messages. Default=TRUE.
#' 
#' @return return GBS.array object with kid hapltoypes in slot \code{gbs_kids}. 
#'   
#'   See \url{https://github.com/yangjl/imputeR/blob/master/vignettes/imputeR-vignette.pdf} for more details.
#'   
#' @examples
#' 
#' set.seed(123457)
#' GBS.array <- sim.array(size.array=50, numloci=10000, hom.error = 0.02, het.error = 0.8,
#'                        rec = 0.25, selfing = 0.5, imiss = 0.5, misscode = 3)
#' #get perfectly phased parents
#' GBS.array <- get_phased(GBS.array, maxchunk=3)
#' #get probability matrices
#' probs <- error_mx(hom.error=0.02, het.error=0.8, imiss=0.5) 
#' phased <- impute_kid(GBS.array, winsize=10, verbose=TRUE)
#' 
#' 
impute_kid <- function(geno, pp, ped, kid_idx=1:10, verbose=TRUE){
    
    #ped <- GBS.array@pedigree
    ped[, 1:3] <- apply(ped[, 1:3], 2, as.character)
    ## run one kid
    ig <- lapply(kid_idx, function(i){
        if(verbose){ message(sprintf("###>>> start to impute kid [ %s ] ...", i )) }
        subgeno_all <- geno[, c("snpid", ped$proid[i])]
        pp1 <- pp[[ped$parent1[i]]]
        pp2 <- pp[[ped$parent2[i]]]
        
        if(sum(pp1$snpid != subgeno_all$snpid) > 0 | sum(pp2$snpid != subgeno_all$snpid) > 0){
            stop("SNPID not in the same order!")
        }
        khap <- one_kid_hap(pp1, pp2, subgeno_all)
        
        ### return the genotype
        return_geno(khap, subgeno_all, pp1, pp2)
    })
    
    if(verbose){ message(sprintf("###>>> Prepare data.frame for output!")) }
    outg <- ig[[1]]
    for(i in 2:length(ig)){
        outg <- merge(outg, ig[[i]], by="snpid")
    }
    return(outg)
}

return_geno <- function(khap, subgeno_all, pp1, pp2){
    
    a <- pp1$hap1 + pp2$hap1
    a[a>2] <- 3
    subgeno_all[,2] <- a
    subgeno_all[subgeno_all$snpid %in% khap$snpid, 2] <- khap$geno
    
    return(subgeno_all)
}

#' @rdname impute_kid
one_kid_hap <- function(pp1, pp2, subgeno_all){
    
    chr <- data.frame()
    for(i in 1:10){
        p1 <- subset(pp1, chr==i)
        p2 <- subset(pp2, chr==i)
        
        cs <- unique(p1$chunk)
        mychunk <- lapply(1:length(cs), function(x){
            print(x)
            hap_in_chunk(p1, p2, c=cs[x], subgeno=subset(subgeno_all, snpid %in% p1$snpid))
        })
        for(j in 1:length(mychunk)){
            chr <- rbind(chr, mychunk[[j]])
        }
    }
    chr$geno <- chr$hap1 + chr$hap2
    return(chr)
}

#' @rdname impute_kid
hap_in_chunk <- function(p1, p2, c, subgeno){
    
    mysnpid <- subset(p1, chunk == c)$snpid
    p1chunk <- subset(p1, snpid %in% mysnpid & hap1 != hap2)
    p2chunk <- subset(p2, snpid %in% mysnpid & hap1 != hap2)
    hetsites <- sort(unique(c(p1chunk$idx, p2chunk$idx)))
    
    ### make sure there is no missing data for the hetsites
    idx1 <- subset(p1, idx %in% hetsites & hap1 ==3)$idx
    idx2 <- subset(p2, idx %in% hetsites & hap1 ==3)$idx
    hetsites <- sort(unique(c(p1chunk$idx, p2chunk$idx)))
    if(length(c(idx1, idx2)) > 0){
        hetsites <- hetsites[!(hetsites %in% c(idx1, idx2))]
    }
    if(length(hetsites) > 0 ){
        
        if(length(unique(p2chunk$chunk)) <= 10){
            hetsnpid <- subset(p1, idx %in% hetsites)$snpid
            ### het sites < window length
            p1_haps <- list(p1[hetsites,]$hap1, p1[hetsites,]$hap2)
            p2_haps <- setup_dad_haps(df=p2[hetsites,], hapcol=4)
            khaps <- which_kid_hap(p1_haps, p2_haps, kidwin=subgeno[subgeno$snpid %in% hetsnpid,2])
            mychunk <- data.frame(hap1=khaps[[1]], hap2=khaps[[2]], snpid=hetsnpid)
            return(mychunk) 
        }else{
            hap_in_large_chunk(p1c=p1[hetsites,], p2c=p2[hetsites,])
        }
        
    }else{
        return(NULL)
    }
}


hap_in_large_chunk <- function(p1c=p1[hetsites,], p2c=p2[hetsites,], subgeno){
    
    cs <- unique(p2c$chunk)
    mychunk <- lapply(1:length(cs), function(x){
        mysnpid <- subset(p2c, chunk == cs[x])$snpid
        p1chunk <- subset(p1c, snpid %in% mysnpid & hap1 != hap2)
        p2chunk <- subset(p2c, snpid %in% mysnpid & hap1 != hap2)
        hetsites <- sort(unique(c(p1chunk$idx, p2chunk$idx)))
        
        ### make sure there is no missing data for the hetsites
        idx1 <- subset(p1, idx %in% hetsites & hap1 ==3)$idx
        idx2 <- subset(p2, idx %in% hetsites & hap1 ==3)$idx
        hetsites <- sort(unique(c(p1chunk$idx, p2chunk$idx)))
        if(length(c(idx1, idx2)) > 0){
            hetsites <- hetsites[!(hetsites %in% c(idx1, idx2))]
        }
        
        if(length(hetsites) > 0 ){
            hetsnpid <- subset(p1c, idx %in% hetsites)$snpid
            ### het sites < window length
            p2_haps <- list(p2c[p2c$snpid %in% hetsnpid,]$hap1, p2c[p2c$snpid %in% hetsnpid,]$hap2)
            p1_haps <- setup_dad_haps(df=p1c[p1c$snpid %in% hetsnpid,], hapcol=4)
            khaps <- which_kid_hap(p1_haps, p2_haps, kidwin=subgeno[subgeno$snpid %in% hetsnpid,2])
            mychunk <- data.frame(hap1=khaps[[1]], hap2=khaps[[2]], snpid=hetsnpid)
            return(mychunk)
        }else{
            return(NULL)
        }
    })
    temchr <- data.frame()
    for(j in 1:length(mychunk)){
        temchr <- rbind(temchr, mychunk[[j]])
    }
    return(temchr)
}

# give this mom haplotype and a kid's diploid genotype over the window and returns maximum prob
# Mendel is takenh care of in the probs[[]] matrix already
#' @rdname impute_kid
#' @param p1_haps Two haplotypes in a chunk of the first parent.
#' @param p2_haps All possible haplotypes in chunks of the 2nd parent.
#' @param kidwin The genotypes of the kid for all heterozygote site of both parents.
which_kid_hap <- function(p1_haps, p2_haps, kidwin){
    genotypes=list()
    #haplotype=unlist(haplotype)
    tab <- data.frame()
    p <- 0
    for(i in 1:length(p1_haps)){
        for(j in 1:length(p2_haps)){
            p <- p +1
            tem <- data.frame(idx1=i, idx2=j, tot=p)
            tab <- rbind(tab, tem)
            genotypes[[p]] <- p1_haps[[i]] + p2_haps[[j]]
        }
    }
    
    geno_probs=as.numeric() #prob of each of three genotypes
    for(geno in 1:p){
        #log(probs[[2]][three_genotypes,kidwin] is the log prob. of kid's obs geno 
        #given the current phased geno and given mom is het. (which is why probs[[2]])
        h1 <- p1_haps[[tab$idx1[geno]]]
        h2 <- p2_haps[[tab$idx2[geno]]]
        geno_probs[geno]=sum( unlist(lapply(1:length(h1), function(zz) 
            log( probs[[h1[zz]+1]][[h2[zz]+1]][genotypes[[geno]][zz]+1, kidwin[zz]+1]))))
    }
    
    maxidx <- which.max(geno_probs)
    return(list(p1_haps[[tab$idx1[maxidx]]], p2_haps[[tab$idx2[maxidx]]]))
}

#' @rdname phase_chunk
#' @param temmom Must be a data.frame(hap1, hap2, chunk, idx)
setup_dad_haps <- function(df, hapcol=4){
    haps1 <- setup_haps(length(unique(df$chunk)))
    haps2 <- lapply(1:length(haps1), function(x) 1-haps1[[x]])
    allhaps <- c(haps1, haps2)
    mom_haps <- lapply(1:length(allhaps), function(x){
        temout <- c()
        k = 1
        for(c in unique(df$chunk)){
            temout <- c(temout, df[df$chunk==c, allhaps[[x]][k]+hapcol])
            k <- k+1
        }
        return(temout)    
    })
    return(mom_haps)
}

