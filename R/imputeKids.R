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
impute_kid <- function(GBS.array, winstart=10, winend=500, stepsize=10, recomb=1.5, verbose=TRUE){
    
    ped <- GBS.array@pedigree
    progeny <- GBS.array@gbs_kids
    phasedp <- GBS.array@gbs_parents
    
    for(k in 1:nrow(ped)){
        kid <- progeny[[ped$kid[k]]]
        p1 <- phasedp[[ped$p1[k]]]
        p2 <- phasedp[[ped$p2[k]]]
        
        if(verbose){ message(sprintf("###>>> start to impute kid [ %s ] ...", k )) }
        win <- optimal_win(p1, p2, kid, winstart, winend, stepsize, recomb)
        if(verbose){ message(sprintf("###>>> optimal win length [ %s ] SNPs, diff1=%s, diff2=%s ...", 
                                     win$win_length, win$bp1, win$bp2)) }
        
        kidgeno <- data.frame()
        kidgeno <- impute_onekid(momphase, kid, win_length=win$win_length, verbose=FALSE)
        progeny[[k]][[3]] <- kidgeno
    }
    return(progeny)
}

#' \code{Maximum likelihood of the inherited haplotypes}
#' 
#' Determining the optimal window size by searching from winstart to winend.
#' And then calculate the maximum likelihood of the inherited haplotypes in the determined window size. 
#'
#' @param GBS.array A GBS.array object with phased parents. 
#' It can be generated from \code{phase_parent} or\code{sim.array} functions.
#' 
#' @param winstart The min window size for the heterozygote sites of the two parents.
#' @param winend The max window size for the heterozygote sites of the two parents.
#' @param stepsize The step size.
#'  
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
optimal_win <- function(p1, p2, kid, winstart, winend, stepsize, recomb){
    bps <- data.frame()
    for(win_length in seq(from=winstart, to=winend, by=stepsize)){
        
        
        
        ###???
        for(c in unique(p1$chunk)){
            ### find the best haps in a chunk
            mychunk <- hap_in_chunk(p1, p2, c, win_length, kid)
            
            ### refine break point
            out <- refine_brk_point(mychunk, win_length, verbose=FALSE, kid)
            #mychunk <- out[[1]]
            tem <- data.frame(win_length=win_length, bp1= tem$bp1+ length(out[[2]][[1]]), 
                              bp2= tem$bp2+ length(out[[2]][[2]]))
        }
        bps <- rbind(bps, tem)
        
    }
    bps$bp1 <- abs(bps$bp1 - recomb)
    bps$bp2 <- abs(bps$bp2 - recomb)
    idx <- which.min(bps$bp1+bps$bp2)
    return(bps[idx,])
}

#' @rdname optimal_win
hap_in_chunk <- function(p1, p2, c, win_length, kid){
    
    p1chunk <- subset(p1, chunk == c & hap1 != hap2)
    p2chunk <- subset(p2[subset(p1, chunk == c)$idx, ], hap1 != hap2)
    hetsites <- sort(unique(c(p1chunk$idx, p2chunk$idx)))
    
    ### het sites < window length
    if(win_length >= length(hetsites)){
        p1_haps <- list(p1[hetsites,]$hap1, p1[hetsites,]$hap2)
        p2_haps <- setup_mom_haps(temmom=p2[hetsites,])
        khaps <- which_kid_hap(p1_haps, p2_haps, kidwin=kid[hetsites])
    }else{ # het sites > window length, use non-overlap window
        for(win in 1:floor(length(hetsites)/win_length) ){
            myidx <- ((win-1)*win_length+1) : (win*win_length)
            myhetsites <- hetsites[myidx]
            p1_haps <- list(p1[myhetsites,]$hap1, p1[myhetsites,]$hap2)
            p2_haps <- setup_mom_haps(temmom=p2[myhetsites,])
            khaps <- which_kid_hap(haplotype, kidwin=kid[myhetsites])
        }
        ##### calculate last window
        myidx <- (length(hetsites)-win_length+1) : length(hetsites)
        myhetsites <- hetsites[myidx]
        p1_haps <- list(p1[myhetsites,]$hap1, p1[myhetsites,]$hap2)
        p2_haps <- setup_mom_haps(temmom=p2[myhetsites,])
        khaps <- which_kid_hap(haplotype, kidwin=kid[myhetsites])   
    }
    mychunk <- data.frame(hap1=khaps[[1]], hap2=khaps[[2]], chunk=c, idx=hetsites)
    return(mychunk)
}

# give this mom haplotype and a kid's diploid genotype over the window and returns maximum prob
# Mendel is takenh care of in the probs[[]] matrix already
#' @rdname optimal_win
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





#' @rdname optimal_win
refine_brk_point <- function(mychunk, win_length, verbose, kid){
    
    bpoints <- get_break_point(mychunk)
    if(verbose){ message(sprintf("###>>> observed number of break points [ %s ] and [ %s ]", 
                                 length(bpoints[[1]]), length(bpoints[[2]]) )) } 
    
    #kidwin <- kid[mychunk$idx]
    mychunk <- move_bp(mychunk, whichhap=1, idx=bpoints[[1]], win_length, kid)
    mychunk <- move_bp(mychunk, whichhap=2, idx=bpoints[[2]], win_length, kid)
    
    bpoints <- get_break_point(mychunk)
    if(verbose){ message(sprintf("###>>> After refining: break points [ %s ] and [ %s ]", 
                                 length(bpoints[[1]]), length(bpoints[[2]]) )) } 
    
    return(list(mychunk, bpoints))
}

#' @rdname optimal_win
get_break_point <- function(mychunk){
    mychunk$r1 <- 0
    # compute the minimum distance to two haplotypes
    mychunk[mychunk$k1!=3,]$r1 <- ifelse(mychunk[mychunk$k1!=3,]$k1 == mychunk[mychunk$k1!=3, ]$hap1, 1, 2)
    mychunk$r2 <- 0
    mychunk[mychunk$k2!=3,]$r2 <- ifelse(mychunk[mychunk$k2!=3,]$k2 == mychunk[mychunk$k2!=3, ]$hap1, 1, 2)
    
    x1 <- factor(paste0(head(mychunk$r1,-1), tail(mychunk$r1,-1)), levels = c('11','12','21','22'))
    tab1 <- table(x1)
    x2 <- factor(paste0(head(mychunk$r2,-1), tail(mychunk$r2,-1)), levels = c('11','12','21','22'))
    tab2 <- table(x2)
    
    #if(sum(tab1[2:3])>2 & sum(tab2[2:3])>2)
    tx <- sum(tab1[2:3])+sum(tab2[2:3])
    idx1 <- sort(c(which(x1=="12"), which(x1=="21")))
    idx2 <- sort(c(which(x2=="12"), which(x2=="21")))
    
    return(list(idx1, idx2))
}



#' \code{Impute kid genotypes}
#'
#' Calculate the maximum likelihood of the inherited haplotypes. 
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
impute_onekid <- function(momphase, kid, win_length, verbose){
    
    res <- data.frame()
    for(c in unique(momphase$chunk)){
        ### find the best haps in a chunk
        mychunk <- hap_in_chunk(momphase, c, win_length, kid)
        #### find the min path of recombinations: it is already minimum path not necessary to run the following
        #mychunk <- min_path(mychunk, verbose)
        
        ### refine break point
        mychunk <- refine_brk_point(mychunk, win_length, verbose, kid)[[1]]
        #mychunk <- out[[1]]
        res <- rbind(res, mychunk)
    }
    return(res)
}






move_bp <- function(mychunk, whichhap, idx, win_length, kid){
    idx <- c(0, idx, nrow(mychunk))
    if(length(idx) > 2){
        for(i in 2:(length(idx)-1)){
            
            ### start and end of the new window
            if((idx[i] - idx[i-1]) > win_length){
                starti <- idx[i] - win_length
            }else{
                starti <- idx[i-1]+1  
            }
            if((idx[i+1] - idx[i]) >= win_length){
                endi <- idx[i] + win_length
            }else{
                endi <- idx[i+1]
            }
            
            ########################################
            haps <- c(mychunk[starti:idx[i], 4+whichhap], 1-mychunk[(idx[i]+1):endi, 4+whichhap])
            hap0 <- lapply(1:(length(haps)-1), function(i) c(haps[1:i], 1-haps[(i+1):length(haps)]))
            hap0[[length(hap0)+1]] <- haps
            hap0[[length(hap0)+1]] <- 1-haps
            fixhap <- mychunk[starti:endi, 7-whichhap]
            genok <- kid[mychunk[starti:endi, ]$idx]
            newhap <- which_fine_hap(hap0, fixhap, genok)
            if(!is.null(newhap)){
                mychunk[starti:endi, 4+whichhap] <- newhap
            }  
        }
    }
    return(mychunk)
}


which_fine_hap <- function(hap0, fixhap, genok){
    #genotype <- hap1+hap2
    geno_probs=as.numeric() #prob of each of three genotypes
    
    for(geno in 1:length(hap0)){
        genotype <- hap0[[geno]] + fixhap
        #log(probs[[2]][three_genotypes,kidwin] is the log prob. of kid's obs geno 
        #given the current phased geno and given mom is het. (which is why probs[[2]])
        geno_probs[geno]=sum( sapply(1:length(genotype), function(zz) log( probs[[2]][genotype[zz]+1, genok[zz]+1])))
    }
    if(length(which(geno_probs==max(geno_probs)))==1){
        return(hap0[[which.max(geno_probs)]])
    }else{
        return(NULL)
    }
}



######################################################




###############################################################
min_path <- function(mychunk, verbose){
    
    mychunk$r1 <- 0
    # compute the minimum distance to two haplotypes
    mychunk[mychunk$k1!=3,]$r1 <- ifelse(mychunk[mychunk$k1!=3,]$k1 == mychunk[mychunk$k1!=3, ]$hap1, 1, 2)
    mychunk$r2 <- 0
    mychunk[mychunk$k2!=3,]$r2 <- ifelse(mychunk[mychunk$k2!=3,]$k2 == mychunk[mychunk$k2!=3, ]$hap1, 1, 2)
    x1 <- factor(paste0(head(mychunk$r1,-1), tail(mychunk$r1,-1)), levels = c('11','12','21','22'))
    tab1 <- table(x1)
    x2 <- factor(paste0(head(mychunk$r2,-1), tail(mychunk$r2,-1)), levels = c('11','12','21','22'))
    tab2 <- table(x2)
    
    #if(sum(tab1[2:3])>2 & sum(tab2[2:3])>2)
    tx <- sum(tab1[2:3])+sum(tab2[2:3])
    idxs <- sort(unique(c(which(x1=="12"), which(x1=="21"), which(x2=="12"), which(x2=="21"))))
    
    if(length(idxs) == 1 ){
        myidx <- (idxs[1]+1):nrow(mychunk)
        out <- compute_transition(mychunk, myidx, tx)
        mychunk <- out[[1]]
        tx <- out[[2]]
    }else if(length(idxs) > 1){
        for(i in 1:(length(idxs)-1)){
            myidx <- (idxs[i]+1):idxs[i+1]
            out <- compute_transition(mychunk, myidx, tx)
            mychunk <- out[[1]]
            tx <- out[[2]]
        }
        myidx <- (idxs[length(idxs)]+1):length(mychunk)
        out <- compute_transition(mychunk, myidx, tx)
        mychunk <- out[[1]]
        tx <- out[[2]]
    }
    
    return(mychunk)    
}
compute_transition <- function(mychunk, myidx, tx){
    mychunk$t1 <- mychunk$r1
    mychunk$t2 <- mychunk$r2
    mychunk$t1[myidx] <- mychunk$r2[myidx]
    mychunk$t2[myidx] <- mychunk$r1[myidx]
    xt1 <- factor(paste0(head(mychunk$t1,-1), tail(mychunk$t1,-1)), levels = c('11','12','21','22'))
    xtab1 <- table(xt1)
    xt2 <- factor(paste0(head(mychunk$t2,-1), tail(mychunk$t2,-1)), levels = c('11','12','21','22'))
    xtab2 <- table(xt2)
    if(sum(xtab1[2:3])+sum(xtab2[2:3]) < tx){
        mychunk$r1 <- mychunk$t1
        mychunk$r2 <- mychunk$t2
        tx2 <- sum(xtab1[2:3])+sum(xtab2[2:3])
        tx <- tx2
    }
    return(list(mychunk[, -9:-10], tx))
}
