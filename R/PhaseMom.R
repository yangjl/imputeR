#'
#' \code{Phasing mom's genotype. } 
#'
#' Phasing mom's genotype using confident genotyping data, for example, WGS data, or 
#' using imputed genotype data.
#'
#' @param obs_mom A vector of mom's GBS data. Coded with 0, 1 and 2, which is the copy of alternative alleles. 
#' Missing data should be coded with 3.
#' @param obs_kids A list of vectors of Kid's GBS data. Coded with 0, 1 and 2, which is the copy of alternative alleles. 
#' Missing data should be coded with 3.
#' @param hom.error Homozygous error rate, default=0.02.
#' @param het.error Heterozygous error rate, default=0.8.
#' @param oddratio The cutoff used for determining the log likelihood ratio of the highest and the 2nd highest genotypes. 
#' The oddratio = NULL means to report the most likely genotype. Default value sets to 0.5.
#' @param returnall For function 'genomom', returnall=TRUE will return with all the information; returnall=FALSE will return
#' mom's genotype only.
#' @return return a data.frame with all the log likelihoods. Using function \code{momgeno} to extract mom's genotype. 
#'   
#'   See \url{https://github.com/yangjl/imputeR/blob/master/vignettes/imputeR-vignette.pdf} for more details.
#' @examples
#' obs_mom <- c(0, 0, 0)
#' obs_kids <- list(c(0, 0, 0), c(0, 0, 0), c(0, 0, 0), c(0, 0, 0), c(0, 0, 0), c(0, 0, 0),
#' c(0, 0, 0), c(0, 0, 0), c(0, 0, 0), c(0, 0, 0), c(0, 0, 0), c(0, 0, 0))
#'
#' impute_parent(obs_mom, obs_kids, hom.error=0.02, het.error=0.8)
#' 
#' parentgeno(geno, oddratio=0.6931472, returnall=FALSE)
#' 
#' 

############################################################################
# Get mom's phase
# should return two haps

### link haplotypes
phasing <- function(estimated_mom, progeny, win_length, verbose=FALSE){
    
    mom_haps <- setup_haps(win_length) 
    if(verbose){ message(sprintf("###>>> start to phase hap chunks ...")) }
    haplist <- phase_mom_chuck(estimated_mom, progeny, win_length, verbose, mom_haps)
    
    if(verbose){ message(sprintf("###>>> start to join hap chunks ...")) } 
    if(length(haplist) > 1){
        out <- joint_mom_chunk(haplist, verbose)
        if(verbose){ message(sprintf("###>>> Reduced chunks from [ %s ] to [ %s ]", length(haplist), length(out))) } 
        return(out)
    }
    else{
        return(haplist)
    }
}
##########################################
joint_mom_chunk <- function(haplist, verbose){
    outhaplist <- list(list())
    outhaplist[[1]] <- haplist[[1]] ### store the extended haps: hap1, hap2 and idx
    i <- 1
    for(chunki in 2:length(haplist)){
        if(verbose){ message(sprintf(">>> join chunks [ %s and %s, total:%s] ...", chunki-1, chunki, length(haplist))) }
        # join two neighbor haplotype chunks
        oldchunk <- haplist[[chunki-1]]
        newchunk <- haplist[[chunki]]
        hapidx <- c(oldchunk[[3]], newchunk[[3]])
        haps <- list(c(oldchunk[[1]], newchunk[[1]]), c(oldchunk[[1]], newchunk[[2]]))
        temhap <- link_haps(momwin=hapidx, progeny, haps, returnhap=FALSE)
        if(!is.null(temhap)){
            outold <- outhaplist[[i]][[1]]
            outoldchunk <- outold[ (length(outold)-length(oldchunk[[1]])+1):length(outold)]
            outnewchunk <- temhap[(length(oldchunk[[1]])+1):length(temhap)]
            
            same <- sum(outoldchunk == temhap[1:length(oldchunk[[1]])])
            #same <- sum(mom_phase1[(length(mom_phase1)-8):length(mom_phase1)] == win_hap[1:length(win_hap)-1])
            outhaplist[[i]][[3]] <- c(outhaplist[[i]][[3]], newchunk[[3]])
            if(same == 0){ #totally opposite phase of last window
                #hap2[length(mom_phase2)+1] <- win_hap[length(win_hap)]
                                
                outhaplist[[i]][[1]] <- c(outhaplist[[i]][[1]], 1-outnewchunk)
                outhaplist[[i]][[2]] <- c(outhaplist[[i]][[2]], outnewchunk)
                
            } else if(same== length(oldchunk[[1]]) ){ #same phase as last window
                #mom_phase1[length(mom_phase1)+1] <- win_hap[length(win_hap)]
                
                outhaplist[[i]][[1]] <- c(outhaplist[[i]][[1]], outnewchunk)
                outhaplist[[i]][[2]] <- c(outhaplist[[i]][[2]], 1-outnewchunk)
            } else{
                stop(">>> Extending error !!!")
            }
        } else {
            i <- i +1
            outhaplist[[i]] <- haplist[[chunki]]
        } 
    }
    return(outhaplist)
}
############################################################################
link_haps <- function(momwin, progeny, haps, returnhap=FALSE){  
    # momwin is list of heterozygous sites, progeny list of kids genotypes, 
    # haps list of possible haps,momphase1 is current phased mom for use in splitting ties
    #### function for running one hap ####
    runoverhaps <- function(myhap){
        #iterate over possible haplotypes <- this is slower because setup_haps makes too many haps
        #get max. prob for each kid, sum over kids
        return(sum( sapply(1:length(progeny), function(z) 
            which_phase(haps[myhap], progeny[[z]][[2]][momwin] ))))
    }
    phase_probs <- sapply(1:(length(haps)), function(a) runoverhaps(a) )
    #if multiple haps tie, return two un-phased haps
    if(length(which(phase_probs==max(phase_probs)))>1){
        return(NULL)
    } else {
        return(haps[[which.max(phase_probs)]])
    }
}
##########################################
phase_mom_chuck <- function(estimated_mom, progeny, win_length, verbose, mom_haps){
    hetsites <- which(estimated_mom==1)
    # gets all possible haplotypes for X hets 
    
    mom_phase1 = mom_phase2 = as.numeric() 
    win_hap = old_hap = nophase = as.numeric() 
    haplist <- list()
    #for(winstart in 1:(length(hetsites)-(win_length-1)))
    winstart <- i <- 1
    
    
    ###### print progress bar
    pb <- txtProgressBar(min = winstart, max = length(hetsites)-(win_length-1), style = 3)
    while(winstart <= length(hetsites)-(win_length-1)){
        if(verbose){ setTxtProgressBar(pb, winstart) } 
        momwin <- hetsites[winstart:(winstart+win_length-1)]
        if(winstart==1){ 
            #arbitrarily assign win_hap to one chromosome initially
            win_hap <- infer_dip(momwin,progeny,haps=mom_haps, returnhap=TRUE)
            mom_phase1=win_hap
            mom_phase2=1-win_hap
            idxstart <- 1
        } else{
            win_hap <- infer_dip(momwin, progeny, haps=mom_haps, returnhap=FALSE)
            ### comparing current hap with old hap except the last bp -JLY
            if(!is.null(win_hap)){
                
                same=sum(mom_phase1[(length(mom_phase1)-win_length+2):length(mom_phase1)]==win_hap[1:length(win_hap)-1])
                
                if(same == 0){ #totally opposite phase of last window
                    mom_phase2[length(mom_phase2)+1] <- win_hap[length(win_hap)]
                    mom_phase1[length(mom_phase1)+1] <- 1-win_hap[length(win_hap)]
                } else if(same==(win_length-1) ){ #same phase as last window
                    mom_phase1[length(mom_phase1)+1] <- win_hap[length(win_hap)]
                    mom_phase2[length(mom_phase2)+1] <- 1-win_hap[length(win_hap)]
                } else{
                    diff1 <- sum(abs(mom_phase1[(length(mom_phase1)-win_length+2):length(mom_phase1)]-win_hap[1:length(win_hap)-1]))
                    diff2 <- sum(abs(mom_phase2[(length(mom_phase1)-win_length+2):length(mom_phase1)]-win_hap[1:length(win_hap)-1]))
                    if(diff1 > diff2){ #momphase1 is less similar to current inferred hap
                        mom_phase2[length(mom_phase2)+1] <- win_hap[length(win_hap)]
                        mom_phase1[length(mom_phase1)+1] <- 1-win_hap[length(win_hap)]
                    } else{ #momphase1 is more similar
                        mom_phase1[length(mom_phase1)+1] <- win_hap[length(win_hap)]
                        mom_phase2[length(mom_phase2)+1] <- 1-win_hap[length(win_hap)]
                    }
                }
            } else {
                ### potential recombination in kids, output previous haps and jump to next non-overlap window -JLY###
                idxend <- winstart + win_length -2
                haplist[[i]] <- list(mom_phase1, mom_phase2, hetsites[idxstart:idxend])
                i <- i +1
                
                ### warning(paste("Likely recombination at position", winstart+1, sep=" "))
                ### if new window is still ambiguous, add 1bp and keep running until find the best hap
                winstart <- winstart + win_length -2
                while(is.null(win_hap)){
                    
                    winstart <- winstart + 1
                    win_hap <- jump_win(winstart, win_length, hetsites, mom_haps)
                    if(is.null(win_hap)){
                        nophase <- c(nophase, hetsites[winstart])
                    }
                }
                idxstart <- winstart
                mom_phase1 <- win_hap
                mom_phase2 <- 1-win_hap
            }
        }
        winstart <- winstart + 1
    }
    close(pb)
    ### return the two haplotypes
    #myh1 <- replace(estimated_mom/2, hetsites, mom_phase1)
    #myh2 <- replace(estimated_mom/2, hetsites, 1-mom_phase1)
    #return(data.frame(h1=myh1, h2=myh2))
    #if(verbose){ message(sprintf(">>> phasing done!")) }
    haplist[[i]] <- list(mom_phase1, mom_phase2, hetsites[idxstart:length(hetsites)])
    ## list: hap1, hap2 and idx; info
    return(haplist)
    #return(list(haplist=haplist, info=list(het=hetsites, nophase=nophase)))
}



############################################################################
jump_win <- function(winstart, win_length, hetsites, mom_haps){
    ### jump to next window
    if(length(hetsites) > (winstart + win_length - 1)){
        momwin <- hetsites[winstart:(winstart + win_length - 1)]
        win_hap <- infer_dip(momwin, progeny, haps=mom_haps, returnhap=FALSE)
    }else{
        momwin <- hetsites[winstart:length(hetsites)]
        mom_haps_tem <- setup_haps(win_length=length(winstart:length(hetsites)))
        win_hap <- infer_dip(momwin, progeny, haps=mom_haps_tem, returnhap=TRUE)
        
    }
    return(win_hap)
}


############################################################################
# Setup all possible haplotypes for window of X heterozgous sites
# This needs to be fixed to remove redundancy. E.g. 010 is the same as 101 and 1010 is same as 0101. 
# I don't think should bias things in the meantime, just be slow.
setup_haps <- function(win_length){
    if(win_length <= 20){
        alist <- lapply(1:win_length, function(a) c(0,1) )
        ### give a combination of all 0,1 into a data.frame
        hapdf <- expand.grid(alist)[1:2^(win_length-1),]
        ### split the data.frame into a list
        return(as.list(as.data.frame(t(hapdf)))) 
    }else{
        stop("!!! Can not handle [win_length > 20] !")
    }   
}
#system.time(tem2 <- setup_haps2(10))
#system.time(tem <- setup_haps(10))


############################################################################
# Infer which phase is mom in a window
infer_dip <- function(momwin, progeny, haps, returnhap=FALSE){  
    # momwin is list of heterozygous sites, progeny list of kids genotypes, 
    # haps list of possible haps,momphase1 is current phased mom for use in splitting ties
    #### function for running one hap ####
    runoverhaps <- function(myhap){
        #iterate over possible haplotypes <- this is slower because setup_haps makes too many haps
        #get max. prob for each kid, sum over kids
        return(sum( sapply(1:length(progeny), function(z) 
            which_phase(haps[myhap],progeny[[z]][[2]][momwin] ))))
    }
    phase_probs <- sapply(1:(length(haps)), function(a) runoverhaps(a) )
    #if multiple haps tie, check each against current phase and return one with smallest distance
    if(length(which(phase_probs==max(phase_probs)))>1){
        if(returnhap){
            return(haps[[sample(which(phase_probs==max(phase_probs)), 1)]])
        } else{
            return(NULL)
        }
    } else if(length(which(phase_probs==max(phase_probs)))==1) {
        return(haps[[which.max(phase_probs)]])
    } else{
        return(NULL)
    }
}
############################################################################

############################################################################
# Find most likely phase of kid at a window, return that probability
# give this mom haplotype and a kid's diploid genotype over the window and returns maximum prob
# Mendel is taken care of in the probs[[]] matrix already 
which_phase <- function(haplotype,kidwin){
    three_genotypes=list()
    haplotype=unlist(haplotype)
    three_genotypes[[1]]=haplotype+haplotype
    three_genotypes[[2]]=haplotype+(1-haplotype)
    three_genotypes[[3]]=(1-haplotype)+(1-haplotype)
    geno_probs=as.numeric() #prob of each of three genotypes
    for(geno in 1:3){
        #log(probs[[2]][three_genotypes,kidwin] is the log prob. of kid's obs geno 
        #given the current phased geno and given mom is het. (which is why probs[[2]])
        geno_probs[geno]=sum( sapply(1:length(haplotype), function(zz) 
            log( probs[[2]][three_genotypes[[geno]][zz]+1,kidwin[zz]+1])))
    }
    ### may introduce error
    #if(length(which(geno_probs==max(geno_probs)))!=1){recover()}
    return(max(geno_probs))
}

###########################################
checkphasing <- function(newmom, sim){
    diff <- tot <- 0
    for(i in 1:length(newmom)){
        truehap <- sim[[1]][newmom[[i]][[3]],]
        esthap <- data.frame(h1=newmom[[i]][[1]], h2=newmom[[i]][[2]])
        tab <- cbind(truehap, esthap)
        idx <- which.max(c(cor(tab$hap1, tab$h1), cor(tab$hap1, tab$h2)))
        
        if(idx == 1){
            a <- nrow(subset(tab, hap1 != h1))
        }else{
            a <- nrow(subset(tab, hap1 != h2))
        }
        
        diff <- diff + a
        tot <- tot + nrow(tab)
    }
    message(sprintf(">>> [ %s ] chunks and [ %s ] error rate", length(newmom), diff/tot))
    return(data.frame(piece=length(newmom), diff=diff, tot=tot))
}

write_mom <- function(newmom){
    momdf <- data.frame()
    for(i in 1:length(newmom)){
        tem <- data.frame(chunk=i, idx=newmom[[i]][[3]], hap1=newmom[[i]][[1]], hap2=newmom[[i]][[2]])
        momdf <- rbind(momdf, tem)
    }
    return(momdf)
}