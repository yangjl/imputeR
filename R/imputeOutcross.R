
##############################
imputingOX <- function(momphase, progeny, winstart, winend, stepsize, expect_recomb=1.5, verbose){
    
    for(k in 1:length(progeny)){
        kid <- progeny[[k]][[2]]
        kidgeno <- data.frame()
        
        if(verbose){ message(sprintf("###>>> start to imputekid [ %s ] with [ %s ] chunks ...", 
                                     k, length(unique(momphase$chunk)) )) }
        win <- optimal_win(momphase, kid, winstart, winend, stepsize, expect_recomb)
        if(verbose){ message(sprintf("###>>> optimal win length [ %s ] SNPs, diff1=%s, diff2=%s ...", win$win_length, win$bp1, win$bp2)) }
        
        kidgeno <- impute_onekid(momphase, kid, win_length=win$win_length, verbose=FALSE)
        progeny[[k]][[3]] <- kidgeno
    }
    return(progeny)
}

######################################################
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
optimal_win <- function(momphase, kid, winstart, winend, stepsize, expect_recomb){
    bps <- data.frame()
    
    for(win_length in seq(from=winstart, to=winend, by=stepsize)){
        
        tem <- data.frame(win_length=winstart, bp1=0, bp2=0)
        for(c in unique(momphase$chunk)){
        ### find the best haps in a chunk
            mychunk <- hap_in_chunk(momphase, c, win_length, kid)
            #### find the min path of recombinations: it is already minimum path not necessary to run the following
            #mychunk <- min_path(mychunk, verbose)
            
            ### refine break point
            out <- refine_brk_point(mychunk, win_length, verbose=FALSE, kid)
            #mychunk <- out[[1]]
            tem <- data.frame(win_length=win_length, bp1= tem$bp1+ length(out[[2]][[1]]), bp2= tem$bp2+ length(out[[2]][[2]]))
        }
        bps <- rbind(bps, tem)
           
    }
    bps$bp1 <- abs(bps$bp1 - expect_recomb)
    bps$bp2 <- abs(bps$bp2 - expect_recomb)
    idx <- which.min(bps$bp1+bps$bp2)
    return(bps[idx,])
}


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
hap_in_chunk <- function(momphase, c, win_length, kid){
    mychunk <- subset(momphase, chunk == c)
    
    #winstart <- i <- 1
    ### kids haplotypes by window
    mychunk$k1 <- 3
    mychunk$k2 <- 3
    if(win_length >= nrow(mychunk)){
        #myidx <- mychunk$idx
        haplotype <- mychunk$hap1
        khaps <- which_kid_hap(haplotype, kidwin=kid[mychunk$idx])
        mychunk <- copy_phase(haplotype, mychunk, khaps, idx=1:nrow(mychunk))
    }else{
        for(win in 1:floor(nrow(mychunk)/win_length) ){
            
            myidx <- ((win-1)*win_length+1) : (win*win_length)
            haplotype <- mychunk$hap1[myidx]
            khaps <- which_kid_hap(haplotype, kidwin=kid[mychunk$idx[myidx]])
            mychunk <- copy_phase(haplotype, mychunk, khaps, idx=myidx)
        }
        ##### calculate last window
        myidx <- (nrow(mychunk)-win_length+1) : nrow(mychunk)
        haplotype <- mychunk$hap1[myidx]
        khaps <- which_kid_hap(haplotype, kidwin=kid[mychunk$idx[myidx]])
        mychunk <- copy_phase(haplotype, mychunk, khaps, idx=myidx)    
    }
    return(mychunk)
}
copy_phase <- function(haplotype, mychunk, khaps, idx){
    if(is.null(khaps)){
        mychunk$k1[idx] <- 3
        mychunk$k2[idx] <- 3
    }else if(khaps == 1){
        mychunk$k1[idx] <- haplotype
        mychunk$k2[idx] <- haplotype
    }else if(khaps == 2){
        mychunk$k1[idx] <- haplotype
        mychunk$k2[idx] <- 1-haplotype
    }else if(khaps == 3){
        mychunk$k1[idx] <- 1-haplotype
        mychunk$k2[idx] <- 1-haplotype
    }else{
        stop("###!!! error! Unexpected haplotype value!")
    }
    return(mychunk)
}

# give this mom haplotype and a kid's diploid genotype over the window and returns maximum prob
# Mendel is takenh care of in the probs[[]] matrix already
#  chunk idx hap1 hap2
#1     1   2    1    0
#2     1  19    1    0
#3     1  20    0    1
#4     1  21    0    1
#5     1  24    1    0
#6     1  28    0    1
which_kid_hap <- function(haplotype, kidwin){
    three_genotypes=list()
    #haplotype=unlist(haplotype)
    three_genotypes[[1]]=haplotype+haplotype
    three_genotypes[[2]]=haplotype+(1-haplotype)
    three_genotypes[[3]]=(1-haplotype)+(1-haplotype)
    geno_probs=as.numeric() #prob of each of three genotypes
    for(geno in 1:3){
        #log(probs[[2]][three_genotypes,kidwin] is the log prob. of kid's obs geno 
        #given the current phased geno and given mom is het. (which is why probs[[2]])
        geno_probs[geno]=sum( sapply(1:length(haplotype), function(zz) log( probs[[2]][three_genotypes[[geno]][zz]+1,kidwin[zz]+1])))
    }
    if(length(which(geno_probs==max(geno_probs)))==1){
        return(which.max(geno_probs))
    }else{
        return(which.max(geno_probs))
    }
}

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
