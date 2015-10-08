### methods for the sim object

### format SimSelfer object to JRI's code
sim2input <- function(sim){
    simp <- sim[[1]] # mom
    simk <- sim[[2]] # kids
    progeny <- list()
    for(i in 1:length(simk)){
        ### list of true and observed kids
        progeny[[i]] <- list(simk[[i]][[2]]$hap1+simk[[i]][[2]]$hap2, simk[[i]][[2]]$obs)
    }
    #p <- frq(progeny)
    ### use the perfect parental genotype
    return(list(simp$geno, progeny))
}

simOX_input <- function(sim, n_phased=5, n_chunk=1){
    ### object of sim: output a list of three, [[1]] data.frame of simulated dad [[2]] list of simulated mom
    ### [[3]] list of simulated kids, [[3]][[n=10]], [[[1]] breakpoints of hap1 and hap2 [[2]] data.frame of kid genotype
    
    #unphased_dad <- sim[[1]]$hap1 + sim[[1]]$hap2 
    true_dad <- data.frame(hap1=sim[[1]]$hap1, hap2=sim[[1]]$hap2, geno=sim[[1]]$hap1 + sim[[1]]$hap2)
    if(n_phased > 0){
        mom <- sim[[2]][1:n_phased]
        for(i in 1:n_phased){
            mom[[i]]$chunk <- sort(sample(c(1:n_chunk), nrow(mom[[i]]), replace=TRUE))
        }
    }else{
        mom <- list()
    }
    if(length(sim[[2]]) > (n_phased+1) ) {
        for(i in (n_phased+1):length(sim[[2]])){
            mom[[i]] <- sim[[2]][[i]]$hap1 + sim[[2]][[i]]$hap2
        }
    }
    
    progeny <- list()
    simk <- sim[[3]]
    for(i in 1:length(simk)){
        ### list of true and observed kids
        progeny[[i]] <- list(simk[[i]][[2]]$hap1+simk[[i]][[2]]$hap2, simk[[i]][[2]]$obs)
    }
    
    #### pedigree info
    ped <- data.frame(kid=1:length(progeny), dad=1, mom=1:length(mom))
    
    message(sprintf("###>>> [[1]]: unphased dad (data.frame)"))
    message(sprintf("###>>> [[2]]: phased and unphased mom [ N=%s+%s ] (list of data.frame + vector)", n_phased, length(sim[[2]])-n_phased ))
    message(sprintf("###>>> [[3]]: outcrossed progeny [ N=%s ] (list of list(real, obs))", length(progeny)))
    message(sprintf("###>>> [[4]]: pedigree (data.frame)", length(progeny)))
    return(list(true_dad, mom, progeny, ped))
}

frq <- function(progeny){
    res <- 0
    for(i in 1:length(progeny)){
        res <- res + progeny[[i]][[2]]
    }
    res <- res/(2*length(progeny))
    res <- replace(res, which(res==0), 0.5/length(progeny))
    return(res)
}

get_sim_kids <- function(sim){
    simk <- sim[[2]]
    progeny <- list()
    for(i in 1:length(simk)){
        progeny[[i]] <- simk[[i]][[2]]
    }
    return(progeny)
} 