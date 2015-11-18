#'
#' \code{A progeny array simulator.} 
#'
#' Simulate a S4 GBS.array object with true_parents, gbs_parents, true_kids, gbs_kids and pedigree.
#' 
#' @param size.array Size of the progeny array.
#' @param numloci Number of loci to simulate.
#' 
#' @param hom.error Homozygous error rate, default=0.02.
#' @param het.error Heterozygous error rate, default=0.8.
#' @param rec Recombination rate, default=0.25.
#' @param selfing Proportion of selfed progeny, default=0.5.
#' @param imiss Individual missing rate, default=0.5.
#' @param misscode Missing code, default=3.
#' 
#' @return Return GBS.array object. 
#' Slot1: true_parents, a list of data.frame(hap1, hap2).
#' Slot2: gbs_parents, a list of genotypes. For example, c(1, 2, 2, 0, 3).
#' Slot3: true_kids, a list of data.frame(hap1, hap2).
#' Slot4: gbs_kids, a list of kid genotypes. For example, c(1, 1, 3, 1, 2).
#' Slot5: pedigree, a data.frame (kid, p1, p2). Note, p1 is the focal parent.
#' 
#'   See \url{https://github.com/yangjl/imputeR/blob/master/vignettes/imputeR-vignette.pdf} for more details.
#'   
#' @examples
#' test <- sim.array(size.array=50, numloci=10)
#' class(test)
#' 
sim.array <- function(size.array, numloci, hom.error=0.02, het.error=0.8, rec=0.25, 
                      selfing=0.5, imiss=0.5, misscode=3){
    ### JRI: http://rpubs.com/rossibarra/self_impute
    #############################
    # 1st, make our focal parent
    #############################
    # make neutral SFS
    sfs <- getsfs()
    # get sample of allele freqs from SFS
    p <- sample(sfs, numloci,replace=TRUE)
    # make focal parent using a data.frame
    sim_focal <- data.frame(hap1=ran.hap(numloci, p), hap2=ran.hap(numloci, p))
    
    ################################################################
    # 2nd, get parent array including selfed and outcrossed parents
    ################################################################
    # get a list of outcrossed parents
    outcrossed <- rbinom(n=1, prob=(1-selfing), size=size.array)
    out_parents <- vector("list", outcrossed)
    out_parents <- lapply(1:outcrossed, function(i) 
        data.frame(hap1=ran.hap(numloci,p), hap2=ran.hap(numloci,p)) )
    # now selfed
    self_parents <- vector("list", size.array - outcrossed)
    self_parents <- lapply(1:(size.array-outcrossed), function(i) sim_focal )
    #combine
    if(outcrossed==0){
        parent_array=self_parents  
    }else if (outcrossed==size.array){
        parent_array=out_parents
    }else{
        parent_array=c(out_parents,self_parents)
    }
    #now we make their diploid genotypes, we add the focal parent on to 
    #the end of the parents array
    parents <- lapply(parent_array, function(q) q[,1]+q[,2])
    parents[[size.array+1]] <- c(sim_focal[,1] + sim_focal[,2])
    #finally, add error to make some crappy gbs_parents
    gbs_parents <- lapply(parents, function(a) add_error(a,hom.error,het.error))
    parent_array[[size.array+1]] <- sim_focal
    ################################################################
    # 3rd, make a progeny array for these parents
    ################################################################
    progeny <- vector("list", size.array)
    #use the kid function!
    #each entry in progeny list has two vectors. [[1]] is true genotype, [[2]] is observed
    progeny <- lapply(1:size.array, function(a) 
        kid(p2=list(parent_array[[a]][,1], parent_array[[a]][,2]), 
            p1=list(sim_focal[,1],sim_focal[,2]), het.error, hom.error, rec, imiss, misscode))
    #now setup observed kids
    obs_kids <- true_kids <- list()
    for(i in 1:size.array){ 
        true_kids[[i]] <- progeny[[i]][[1]]
        obs_kids[[i]] <- progeny[[i]][[2]] 
    }
    
    ################################################################
    # 4th, now we impute the focal parent
    ################################################################
    #which parent is our focal one? Here we set to end of parents array for ease
    obs_parent=size.array+1 #focal parent
    #which parents are the other parent of each offspring. These are in order since we simulated them that way.
    #other_parents=c(1:outcrossed,rep(obs_parent,size.array-outcrossed)) #list of other parents
    
    ped <- data.frame(kid=1:size.array, p1= size.array + 1, 
                      p2= c(1:outcrossed,rep(obs_parent,size.array-outcrossed) ))
    
    obj <- new("GBS.array",
               true_parents = parent_array, # list of data.frame(hap1, hap2)
               gbs_parents = gbs_parents,
               true_kids = true_kids,
               gbs_kids = obs_kids,
               pedigree = ped
               )
    return(obj)    
}

#' @rdname sim.array
#' @param p1 first parent (a list)
#' @param p2 second parent (a list)
#' @return a list of true [[1]] and observed [[2]] kid
#'
kid <- function(p1, p2, het.error, hom.error, rec, imiss, misscode){
    if(rec==0){
        k1=p1[[rbinom(1,1,.5)+1]]
        k2=p2[[rbinom(1,1,.5)+1]]
    } else{
        k1=copy.mom(p1,rec) # list
        k2=copy.mom(p2,rec)
    }
    true_kid=k1 + k2
    #return(list(true_kid,obs_kid))
    
    obs_kid <- add_error(true_kid, hom.error, het.error)
    if(imiss > 0){ #don't think the iff statement is necessary
        idx <- missing.idx(length(true_kid), imiss)
        obs_kid <- replace(obs_kid, idx, misscode)
    }
    
    return(list(true_kid, obs_kid))
}

#' @rdname sim.array  
#' @param numloci number of loci.
#' @param p a vector of allele freq with length numloci.
#' @return a vector of genotype,
#' @examples
#' ran.hap(3, c(0.5, 0.8, 0.1))
#'
ran.hap <- function(numloci,p){
    sapply(1:numloci,function(x) rbinom(1,1,p[x]))
}

#' @rdname sim.array  
#' @param diploid a vector of diploid genotype, for example, c(0, 1, 1, 0, 2, 2).
#' @param hom.error homozygous error rate.
#' @param het.error heterozygous error rate.
#' @return a vector of diploid genotype.
#' @examples
#' add_error(c(0, 1, 1, 0, 2, 2), 0.02, 0.8)
#'
add_error <- function(diploid,hom.error,het.error){
    hets_with_error=sample(which(diploid==1),round(het.error*length(which(diploid==1))))
    hom0_with_error=sample(which(diploid==0),round(hom.error*length(which(diploid==0))))
    hom2_with_error=sample(which(diploid==2),round(hom.error*length(which(diploid==2))))
    diploid=replace(diploid,hets_with_error,sample(c(0,2),length(hets_with_error),replace=T)  )
    ### error rate from (hom => het) == (hom1 => hom2)
    diploid=replace(diploid,hom0_with_error,sample(c(1,2),length(hom0_with_error),replace=T)  )
    diploid=replace(diploid,hom2_with_error,sample(c(1,0),length(hom2_with_error),replace=T)  )
    return(diploid)
}

#' @rdname sim.array
#' @param mom a list of mom's haplotype.
#' @param co_mean mean number of crossovers per chromosome.
#' @return a vector of diploid genotype.
#' @examples
#' copy.mom(list(c(0, 1, 1, 0, 1, 1)), 0.25)
copy.mom <- function(mom, co_mean){ 
    co=rpois(1,co_mean) #crossovers
    numloci=length(mom[[1]])
    recp=unique(c(1,sort(round(runif(co, min=2, max=numloci-1))), numloci+1)) #position   
    chrom=rbinom(1,1,.5)+1
    kpiece=as.numeric()
    hap <- c()
    for(r in 1:(length(recp)-1)){
        kpiece=c(kpiece,mom[[chrom]][recp[r]:(recp[r+1]-1)]) #copy 1->rec from mom
        hap <- c(hap, chrom)
        chrom=ifelse(chrom==1,2,1)  
    }     
    return(kpiece)
}

#' @rdname sim.array
missing.idx <- function(nloci, imiss){
    #hist(rbeta(10000, 2, 2))
    if(imiss >= 1){
        m <- rbeta(1, 2, 2)
    }else{
        m <- imiss
    }
    
    ml <- sort(sample(1:nloci, size=round(m*nloci)))
    return(ml)
}

get_true_GBS <- function(GBS.array){
    
    ped <- GBS.array@pedigree
    if(length(unique(ped$p1)) != 1){
        stop("### more than one focal parent!!!")
    }
    
    for(pidx in 1:length(GBS.array@true_parents)){
        true_p <- GBS.array@true_parents[[pidx]]
        GBS.array@gbs_parents[[pidx]] <- true_p$hap1 + true_p$hap2 
    }
    return(GBS.array)
}
