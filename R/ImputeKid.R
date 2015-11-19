#' \code{Phasing and Imputing Kid Genotypes}
#'
#' Calculate the maximum likelihood of the inherited haplotypes. And then find the minimum path of recombinations.
#' At the potential recombination sites, find the maximum likelihood of the breaking points.
#'
#' @param GBS.array A GBS.array object with phased parents. 
#' It can be generated from \code{phase_parent} or\code{sim.array} functions.
#' @param win_length Window length used for halotype phasing. Default=10. 
#' win_length > 20 will dramatically increase computational burden. 
#' @param join_length The length of each neighboring chunks used to connect them into a longer one.
#' @param verbose Writing verbose messages. Default=TRUE.
#' 
#' @return return a data.frame with all the log likelihoods. Using function \code{momgeno} to extract mom's genotype. 
#'   
#'   See \url{https://github.com/yangjl/imputeR/blob/master/vignettes/imputeR-vignette.pdf} for more details.
#'   
#' @examples
#' 
#' GBS.array <- sim.array(size.array=50, numloci=1000, imiss=0.2, selfing=0.5)
#' 
#' 
impute_kid <- function(GBS.array, winsize=10, verbose=TRUE){
    
    for(k in 1:nrow(ped)){
        if(verbose){ message(sprintf("###>>> start to impute kid [ %s ] ...", k )) }
        
        kidhap <- chr_hap(GBS.array, k, winsize)
        
        ped <- GBS.array@pedigree
        GBS.array@gbs_kids[[ped$kid[k]]] <- kidhap
    }
    return(GBS,array)
}

#' @rdname impute_kid
chr_hap <- function(GBS.array, k, winsize){
    ped <- GBS.array@pedigree
    phasedp <- GBS.array@gbs_parents
    kid <- GBS.array@gbs_kids[[ped$kid[k]]]
    p1 <- phasedp[[ped$p1[k]]]
    p2 <- phasedp[[ped$p2[k]]]
    
    chr <- data.frame()
    for(c in unique(p1$chunk)){
        ### find the best haps in a chunk
        mychunk <- hap_in_chunk(p1, p2, c, winsize, kid)
        chr <- rbind(chr, mychunk)
    }
    
    # haps for both homo and heter sites
    hap1 <- p1$hap1
    hap1[chr$idx] <- chr$hap1
    hap2 <- p2$hap2
    hap2[chr$idx] <- chr$hap2
    
    return(data.frame(hap1=hap1, hap2=hap2, chunk1=p1$chunk, chunk2=p2$chunk))
}

#' @rdname impute_kid
hap_in_chunk <- function(p1, p2, c, winsize, kid){
    
    p1chunk <- subset(p1, chunk == c & hap1 != hap2)
    p2chunk <- subset(p2[subset(p1, chunk == c)$idx, ], hap1 != hap2)
    hetsites <- sort(unique(c(p1chunk$idx, p2chunk$idx)))
    
    ### het sites < window length
    if(winsize >= length(hetsites)){
        p1_haps <- list(p1[hetsites,]$hap1, p1[hetsites,]$hap2)
        p2_haps <- setup_mom_haps(temmom=p2[hetsites,])
        khaps <- which_kid_hap(p1_haps, p2_haps, kidwin=kid[hetsites])
    }else{ # het sites > window length, use non-overlap window
        khaps <- list(p1chunk$p1[hetsites], p1chunk$p1[hetsites]) #place holder for the kid haps
        for(win in 1:floor(length(hetsites)/winsize) ){
            myidx <- ((win-1)*winsize+1) : (win*winsize)
            myhetsites <- hetsites[myidx]
            p1_haps <- list(p1[myhetsites,]$hap1, p1[myhetsites,]$hap2)
            p2_haps <- setup_mom_haps(temmom=p2[myhetsites,])
            tem_haps <- which_kid_hap(p1_haps, p2_haps, kidwin=kid[myhetsites])
            
            khaps[[1]][myidx] <- tem_haps[[1]]
            khaps[[2]][myidx] <- tem_haps[[2]]
        }
        ##### calculate last window
        myidx <- (length(hetsites)-winsize+1) : length(hetsites)
        myhetsites <- hetsites[myidx]
        p1_haps <- list(p1[myhetsites,]$hap1, p1[myhetsites,]$hap2)
        p2_haps <- setup_mom_haps(temmom=p2[myhetsites,])
        tem_haps <- which_kid_hap(p1_haps, p2_haps, kidwin=kid[myhetsites])
        
        khaps[[1]][myidx] <- tem_haps[[1]]
        khaps[[2]][myidx] <- tem_haps[[2]]
    }
    mychunk <- data.frame(hap1=khaps[[1]], hap2=khaps[[2]], chunk=c, idx=hetsites)
    return(mychunk)
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
        geno_probs[geno]=sum( sapply(1:length(h1), function(zz) 
            log( probs[[h1[zz]+1]][[h2[zz]+1]][genotypes[[geno]][zz]+1, kidwin[zz]+1])))
    }
    
    maxidx <- which.max(geno_probs)
    return(list(p1_haps[[tab$idx1[maxidx]]], p2_haps[[tab$idx2[maxidx]]]))
}




