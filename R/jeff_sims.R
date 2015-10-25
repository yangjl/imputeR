    
    ### size.array: number of kids
    ### imissing => individual missing rate, 
    ### imiss > 1 will be sampled from a Beta(2,2) distribution (1- U-shaped)
    ### rec: recombination rate 

sim_errors=as.numeric()
for(sim in 1:10){
    misscode = 3
    numloci=500
    hom.error=0.02
    het.error=0.8
    imiss=0.1
    selfing=0.5
    size.array=50
    rec=0.25
    # make mom
    sfs <- getsfs()
    #get freqs for all loci

    p <- sample(sfs, numloci, replace=TRUE) 
    
    ### make focal_parent using a data.frame
    sim_focal <- data.frame(hap1=ran.hap(numloci,p), hap2=ran.hap(numloci,p))
    
    ### start with outcrossed
    outcrossed=rbinom(n=1,prob=selfing,size=size.array)
    out_parents <- vector("list", outcrossed)
    out_parents <- lapply(1:outcrossed, function(i)
            data.frame(hap1=ran.hap(numloci,p), hap2=ran.hap(numloci,p)) 
        )
    #now selfed
    self_parents <- vector("list", size.array-outcrossed)
    self_parents <- lapply(1:(size.array-outcrossed), function(i)
        sim_focal 
    )
    #combine
    parent_array=c(out_parents,self_parents)

    progeny <- vector("list", size.array)
    progeny <- lapply(1:size.array, function(a) 
        kid(p2=list(parent_array[[a]][,1], parent_array[[a]][,2]), p1=list(sim_focal[,1],sim_focal[,2]), 
            het.error, hom.error, rec, imiss, misscode))
    #each entry in progeny list has two vectors. [[1]] is true genotype, [[2]] is observed

    #now setup
    obs_kids=list()
    for(i in 1:size.array){ obs_kids[[i]]=progeny[[i]][[2]] }
    
    parents<-lapply(parent_array, function(q) q[,1]+q[,2]  )
    parents[[size.array+1]]=c(sim_focal[,1]+sim_focal[,2])
    gbs_parents=lapply(parents, function(a) add_error(a,hom.error,het.error))
    
    obs_parent=size.array+1
    other_parents=c(1:outcrossed,rep(obs_parent,size.array-outcrossed))
    inferred_geno_likes=impute_parent(gbs_parents, obs_parent, other_parents, obs_kids, hom.error=0.02, het.error=0.8,p)
    bob=parentgeno(inferred_geno_likes, oddratio=0.6931472, returnall=TRUE)
    sim_errors[sim]=sum(abs(bob$gmax-parents[[obs_parent]]))    
}
    
    
    
#' @rdname SimSelfer Create random haplotype with sfs
ran.hap <- function(numloci,p){
    sapply(1:numloci,function(x) rbinom(1,1,p[x]))
}

#' @rdname SimSelfer 
# Add error to diploid
add_error<-function(diploid,hom.error,het.error){
    hets_with_error=sample(which(diploid==1),round(het.error*length(which(diploid==1))))
    hom0_with_error=sample(which(diploid==0),round(hom.error*length(which(diploid==0))))
    hom2_with_error=sample(which(diploid==2),round(hom.error*length(which(diploid==2))))
    diploid=replace(diploid,hets_with_error,sample(c(0,2),length(hets_with_error),replace=T)  )
    ### error rate from (hom => het) == (hom1 => hom2)
    diploid=replace(diploid,hom0_with_error,sample(c(1,2),length(hom0_with_error),replace=T)  )
    diploid=replace(diploid,hom2_with_error,sample(c(1,0),length(hom2_with_error),replace=T)  )
    return(diploid)
}

#' @rdname SimSelfer 
#Copy mom to kids with recombination
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

#' @rdname SimSelfer 
# add missing
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

# Returns a list of true [[1]] and observed [[2]] kid

#p2=list(parent_array[[a]][,1], parent_array[[a]][,2])
#p1=list(sim_focal[,1],sim_focal[,2])

kid <- function(p1, p2, het.error, hom.error, rec=1.5, imiss=0, misscode=3){
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
    if(imiss > 0){
        idx <- missing.idx(length(true_kid), imiss)
        obs_kid <- replace(obs_kid, idx, misscode)
    }
    
    return(list(true_kid, obs_kid))
}