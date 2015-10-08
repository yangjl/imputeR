### JRI: http://rpubs.com/rossibarra/self_impute
SimOXer <- function(size.array=10, het.error=0.7, hom.error=0.002, numloci=1000, rec=1.5, imiss=0.3, misscode = 3){
    
    ### size.array: number of kids
    ### imissing => individual missing rate, 
    ### imiss > 1 will be sampled from a Beta(2,2) distribution (1- U-shaped)
    ### rec: recombination rate 
    
    # make mom
    sfs <- getsfs()
    #get freqs for all loci
    p <- sample(sfs, numloci, replace=TRUE) 
    
    ### make dad using a data.frame
    sim_dad <- data.frame(hap1=ran.hap(numloci,p), hap2=ran.hap(numloci,p))
    
    ### make an array of mom
    mom_array <- vector("list", size.array)
    mom_array <- lapply(1:size.array, function(i)
        data.frame(hap1=ran.hap(numloci,p), hap2=ran.hap(numloci,p)))
    
    # make selfed progeny array
    progeny <- vector("list", size.array)
    progeny <- lapply(1:size.array, function(a) 
        kid(mom=list(mom_array[[a]][,1], mom_array[[a]][,2]), dad=list(sim_dad[,1],sim_dad[,2]), 
            het.error, hom.error, rec, imiss, misscode))
    #progeny <- replicate(size.array, kid(true_mom,true_mom, het.error, hom.error, recombination=TRUE))
    return(list(sim_dad, mom_array, progeny))
    ### output a list of three, [[1]] data.frame of simulated dad [[2]] list of simulated mom
    ### [[3]] list of simulated kids, [[3]][[n=10]], [[[1]] breakpoints of hap1 and hap2 [[2]] data.frame of kid genotype
}

####
SimSelfer <- function(size.array=20, het.error=0.7, hom.error=0.002, numloci=1000, rec=1.5, imiss=0.3){
    
    ### Simulate and Test
    ### imissing => individual missing rate, 
    ### imiss > 1 will be sampled from a Beta(2,2) distribution (1- U-shaped)
    ### rec: recombination rate 
    
    misscode = 3
    # make mom
    sfs <- getsfs()
    p=sample(sfs,numloci) #get freqs for all loci
    a1=ran.hap(numloci,p) #make haplotypes
    a2=ran.hap(numloci,p)
    
    true_mom=list(a1,a2) #phased 
    #obs_mom=add_error(a1+a2,hom.error,het.error) #convert to diploid genotype
    #if(imiss > 0){
    #    idxmom <- missing.idx(numloci, imiss)
    #    obs_mom <- replace(obs_mom, idxmom, misscode)
    #}
    simp <- data.frame(hap1=a1, hap2=a2, geno=a1+a2, obs=a1+a2)
    
    
    # make selfed progeny array
    progeny <- vector("list",size.array)
    progeny <- lapply(1:size.array, function(a) 
        kid(true_mom, true_mom, het.error, hom.error, rec=rec, imiss=imiss, misscode=misscode))
    #progeny <- replicate(size.array, kid(true_mom,true_mom, het.error, hom.error, recombination=TRUE))
    return(list(simp, progeny))
    
}

### Create random haplotype with sfs
ran.hap <- function(numloci,p){
    sapply(1:numloci,function(x) rbinom(1,1,p[x]))
}
### setup the neutral SFS
getsfs <- function(){
    x=1:99/100 #0.01 bins of freq.
    freq=1/x
    sfs=as.numeric(); 
    for(i in 1:99){sfs=c(sfs,rep(x[i],100*freq[i]))}
    return(sfs)
}
### Add error to diploid
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

# Copy mom to kids with recombination
copy.mom <- function(mom, co_mean){ 
    co=rpois(1,co_mean) #crossovers
    numloci=length(mom[[1]])
    recp=c(1,sort(round(runif(co, min=2, max=numloci-1))), numloci+1) #position   
    chrom=rbinom(1,1,.5)+1
    kpiece=as.numeric()
    hap <- c()
    for(r in 1:(length(recp)-1)){
        kpiece=c(kpiece,mom[[chrom]][recp[r]:(recp[r+1]-1)]) #copy 1->rec from mom
        hap <- c(hap, chrom)
        chrom=ifelse(chrom==1,2,1)  
    }     
    return(list(kpiece, data.frame(hap=hap, start=recp[-length(recp)], end=recp[-1]) ))
}

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

############################################################################
# Make a kid
#Returns list of true [[1]] and observed [[2]] kid
kid <- function(mom, dad, het.error, hom.error, rec=1.5, imiss=0.3, misscode=3){
    if(rec==0){
        k1=mom[[rbinom(1,1,.5)+1]]
        k2=dad[[rbinom(1,1,.5)+1]]
    } else{
        k1=copy.mom(mom,rec) # list
        k2=copy.mom(dad,rec)
    }
    true_kid=k1[[1]] + k2[[1]]
    #return(list(true_kid,obs_kid))
    
    obs_kid <- add_error(true_kid, hom.error, het.error)
    if(imiss > 0){
        idx <- missing.idx(length(true_kid), imiss)
        obs_kid <- replace(obs_kid, idx, misscode)
    }
    
    info <- list(k1[[2]], k2[[2]])
    simk <- data.frame(hap1=k1[[1]], hap2=k2[[1]], obs= obs_kid )
    return(list(info, simk))
}