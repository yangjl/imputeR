#'
#' \code{Create GBS.array objects} 
#'
#' Create GBS.array objects from Geno4imputeR objects. It will compute population allele frq based on families
#' of selfed kids. Original Geno4imputeR will be updated and returned. A set of GBS.array objects will be output
#' to the directory by chromosomes.
#' 
#' @param geno a data.frame of gentoype. coded with 0 (major), 1 (het), 2 (minor) and 3 (missing).
#' @param ped a data.frame of pedigree information, must include three columns: proid (progeny id),
#' parent1 (the first parent id), parent2 (the 2nd parent id). 
#' @param pargeno A data.frame with two columns: "parentid" and "true_p". 
#' The 2nd column indicates whether the true parental genotype is available.
#' @param outdir path of the output dir.
#' @param bychr Wether to split the data by 10 chromsomes.
#' @param snpinfo Wether to calculate snpinfo (this data.frame can be obtained by using function of get_snpinfo).
#' @param size_cutoff Minimum family size required for calculate snpinfo, default=30.
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
create_array <- function(geno, ped, pargeno, pp=NULL, pinfo=NULL, 
                         outdir="largedata/obs", bychr=FALSE, bychunk=NULL){
    
    ### check genotype and pedigree data
    if( sum(!unique(c(ped$proid, ped$parent1, ped$parent2)) %in% names(geno)) > 0 ){
        stop("###>>> Some plant ID could not be found in the genotype file !")
    }else{
        message(sprintf("###>>> Loaded [ %s ] biallelic loci for [ %s ] plants", nrow(geno), nrow(ped)))
    }
    
    #### get pedigree info
    if(is.null(pinfo)){
        pinfo <- pedinfo(ped)
        message(sprintf("###>>> Calculating pedigree info ..."))
    }else{
        message(sprintf("###>>> Loaded pedigree info ..."))
    }
    #### get snpinfo
    #if(is.null(snpinfo)){
    #    snpinfo <- get_snpinfo(geno, ped, self_cutoff)
    #    message(sprintf("###>>> Calculating SNP info ..."))
    #}else{
    #    message(sprintf("###>>> Loaded SNP info ..."))
    #}
    
    message(sprintf("###>>> Preparing GBS.array objects into [ %s ] chr and [ %s ] chunks", 
                    10*sum(bychr), sum(bychunk)))
    #geno <- merge(snpinfo, geno[, -2:-3], by="snpid")
    geno$major <- as.numeric(as.character(gsub("^S|_.*$", "", geno$snpid)))
    geno$minor <- as.numeric(as.character(gsub("^S.*_", "", geno$snpid)))
    names(geno)[2:3] <- c("chr", "pos")
    for(i in 1:nrow(pinfo)){
        ### get pedigree and idx for the p1 and p2
        focalp <- as.character(pinfo$founder[i])
        myped <- ped_focal(ped, pargeno, focalp)
        
        message(sprintf("###>>> Preparing for the [ %sth ] focal parent [ %s ]", i, focalp)) 
        message(sprintf("###>>> It has [ %s selfed ] + [ %s outcrossed/ %s know ] kids ... ", 
                        nrow(subset(myped, p1==p2)), nrow(subset(myped, p1 != p2)), sum(myped$true_p) ))

        ### get snp matrix of each chr
        if(bychr){
            for(chrj in 1:10){
                subgeno <- subset(geno, chr == chrj)
                #subgeno <- subgeno[order(subgeno$pos),]
                
                ### divided into chunks by hetsites
                if(!is.null(bychunk)){
                    
                    cks <- divid_chunks(subgeno, bychunk, focalp)
                    message(sprintf("###>>> Divided Chr [ %s ] into [ %s ] chunks!", chrj, length(cks)))
                    
                    for(k in length(cks)){
                        obj <- get_array_item(subgeno, myped, focalp, pp)
                        outfile <- paste0(outdir, "/", focalp,"_chr", chrj, "_chunk", k, ".RData"  )
                        save(file=outfile, list="obj")
                    }
                    
                }else{
                    obj <- get_array_item(subgeno, myped, focalp, pp)
                    outfile <- paste0(outdir, "/", focalp,"_chr", chrj, ".RData"  )
                    save(file=outfile, list="obj")
                }
            }
        }else{
            obj <- get_array_item(subgeno=geno, myped, focalp, pp)
            outfile <- paste0(outdir, "/", focalp,"_chrall", ".RData"  )
            save(file=outfile, list="obj")
        }
    }
    message("###>>> DONE! <<<###")
}


divid_chunks <- function(subgeno, bychunk, focalp){
    
    idx <- which(subgeno[, focalp] ==1)
    totc <- floor(length(idx)/bychunk)
    
    chgeno <- list()
    if(totc <= 1){
        chgeno[[1]] <- subgeno
    }else{
        for(ci in 1:totc){
            if(ci == 1){
                btidx <- idx[bychunk] 
                chgeno[[ci]] <- subgeno[1:btidx,]
                
            }else if(totc > 1 & ci != totc){
                tpidx <- btidx +1
                btidx <- idx[bychunk*ci]
                chgeno[[ci]] <- subgeno[tpidx:btidx,]
            }else if(totc >1 & ci == totc){
                tpidx <- btidx +1
                chgeno[[ci]] <- subgeno[tpidx:nrow(subgeno),]
            }else{
                stop("!!! divid chunk error!!!")
            }
        }
    }
    return(chgeno)
    
}



#' @rdname create_array
get_snpinfo <- function(geno, ped, self_cutoff){

    message("###>>> Calculating allele frq and missing using all data ... ", appendLF=FALSE)
    snpinfo <- geno[,1:3]
    snpinfo$chr <- as.numeric(as.character(gsub("S|_.*", "", geno$snpid)))
    snpinfo$pos <- as.numeric(as.character(gsub(".*_", "", geno$snpid)))
    snpinfo$miss <- est_missing(geno[,-1:-3])
    snpinfo$totmaf <- est_maf(geno[,-1:-3])
    message("done.")
    
    pinfo <- pedinfo(ped)
    pinfo <- pinfo[order(pinfo$tot, decreasing=TRUE),]
    
    message("###>>> Calculating pop frq and missing using selfed family with [ > %s ] kids ... ", self_cutoff)
    snpinfo$frq <- est_pop_maf(geno, pinfo, ped, self_cutoff)
    return(snpinfo)
}

#' @rdname create_array
pedinfo <- function(ped){
    ped[,1:3] <- apply(ped[,1:3], 2, as.character)
    ped$parent1 <- as.character(ped$parent1)
    ped$parent2 <- as.character(ped$parent2)
    
    sf <- subset(ped, parent1 == parent2)
    ox <- subset(ped, parent1 != parent2)
    
    pinfo <- data.frame(table(sf$parent1))
    names(pinfo) <- c("founder", "nselfer")
    
    oxinfo <- data.frame(table(c(ox$parent1, ox$parent2)))
    if(nrow(oxinfo) > 0){
        names(oxinfo) <- c("founder", "nox")
        pinfo <- merge(pinfo, oxinfo, by="founder", all=TRUE)
        #pinfo$sid <- gsub("..:.*", "", pinfo$founder)
        #pinfo$sid <- gsub("_m", "", pinfo$sid)
        
        pinfo[is.na(pinfo)] <- 0
    }else{
        pinfo$nox <- 0
    }

    pinfo$tot <- pinfo$nselfer + pinfo$nox
    message(sprintf("###>>> Detected [ %s ] parents with [ %s/%s ] kids/haps",
                    nrow(pinfo), nrow(ped), sum(pinfo$nselfer, na.rm=T)*2 + sum(pinfo$nox) ))
    return(pinfo[order(pinfo$tot, decreasing=TRUE),])
}

#' @rdname create_array
est_pop_maf <- function(geno, pinfo, ped, self_cutoff){
    self <- subset(pinfo, nselfer > self_cutoff)
    subkids <- subset(ped, parent1 == parent2 & parent1 == as.character(self$founder)[1] )
    
    frq <- geno[, 1:3]
    for(i in 1:nrow(self)){
        kid_id <- subset(ped, parent1 == parent2 & parent1 == as.character(self$founder)[i] )$proid
        kidgeno <- geno[,as.character(kid_id)]
        frq$pop <- est_maf(pgeno=kidgeno)
        names(frq)[ncol(frq)] <- paste0("p", i)
    }
    frq$popfrq <- apply(frq[, -1:-3], 1, function(x) mean(x, na.rm=T))
    minp <- 1/(2*nrow(self))
    if(nrow(frq[frq$popfrq < minp, ]) > 0) frq[frq$popfrq < minp, ]$popfrq <- minp
    return(frq$popfrq)
}

#' @rdname create_array
est_maf <- function(pgeno){
    freq <- apply(pgeno, 1, function(x){
        c0 <- sum(x == 0)
        c1 <- sum(x == 1)
        c2 <- sum(x == 2)
        return((2*c2 + c1)/(2*(c0 + c1 + c2)))
    })
    return(freq)
}

#' @rdname create_array
est_missing <- function(pgeno){
    freq <- apply(pgeno, 1, function(x){
        c3 <- sum(x == 3)
        return(c3/length(x))
    })
    return(freq)
}
#' @rdname create_array
ped_focal <- function(ped, pargeno, focalp){
    myped <- subset(ped, parent1 == focalp | parent2 == focalp)[, 1:3]
    myped$p1 <- 1
    df <- data.frame(parent2=unique(c(focalp, as.character(myped$parent1), as.character(myped$parent2))), p2=1)
    df$p2 <- 1:nrow(df)
    for(p in 1:nrow(myped)){
        if(myped$parent1[p] != focalp){
            myped$parent2[p] <- myped$parent1[p]
            myped$parent1[p] <- focalp
        }  
    }
    myped <- merge(myped, df, by="parent2")
    myped <- myped[, c("proid", "parent1", "parent2", "p1", "p2")]
    
    myped <- merge(myped, pargeno, by.x="parent2", by.y="parentid")
    myped <- myped[, c("proid", "parent1", "parent2", "p1", "p2", "true_p")]
    myped <- myped[order(myped$p2),]
    return(myped)
}

#' @rdname create_array
get_array_item <- function(subgeno, myped, focalp, pp){
   
    gbsp <- gbsk <- list()
    gbsp[[1]] <- as.vector(subgeno[, focalp])
    for(k in 1:nrow(myped)){
        pk <- myped$p2[k]
        if(!is.null(pp) & sum(myped$parent2[k] %in% names(pp)) > 0){
            subpp <- pp[[myped$parent2[k]]]
            subpp <- subpp[subpp$snpid %in% subgeno$snpid,]
            if(sum(names(subpp$snpid) != names(subgeno$snpid)) > 0){
                stop("!!! snpid of phased parents != subgneo !!!")
            }
            gbsp[[pk]] <- subpp
        }else{
            gbsp[[pk]] <- as.vector(subgeno[, myped$parent2[k]])
        }
        gbsk[[k]] <- as.vector(subgeno[, myped$proid[k]])
    }
    myped$kid <- 1:nrow(myped)
    obj <- new("GBS.array",
               #true_parents = list(), # list of data.frame(hap1, hap2)
               gbs_parents = gbsp,
               #true_kids = true_kids,
               gbs_kids = gbsk,
               pedigree = myped,
               snpinfo = subgeno[,1:3])
    return(obj)
}

#'
#' \code{Update GBS.array object} 
#'
#' Update GBS.array object slot gbs_parents.
#' 
#' @param GBS.array input object.
#' @param ipdat a data.frame of imputed gbs parents.
#' 
#' @return Return updated GBS.array object. 
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
#' update_gbs_parents(GBS.array, ipdat)
#' 
update_gbs_parents <- function(GBS.array, ipdat){
    
    ped <- GBS.array@pedigree
    tab <- data.frame(pid=c(ped$parent1, ped$parent2), pidx=c(ped$p1, ped$p2))
    tab <- tab[!duplicated(tab$pidx),]
    tab <- tab[order(tab$pidx),]
    a1 <- nrow(tab)
    a2 <- length(GBS.array@gbs_parents)
    if(a1 != a2) stop("ERROR, unique parents number in pedigree not eq to length of gbs_parents!")
    
    for(idx in 1:nrow(tab)){
        GBS.array@gbs_parents[[idx]] <- ipdat[, tab$pid[idx]]
    }
    return(GBS.array)
}


