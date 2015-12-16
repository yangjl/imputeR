# simulaltion 
true_other_parents <- function(GBS.array){
    for(i in 1:(length(GBS.array@gbs_parents)-1)){
        GBS.array@gbs_parents[[i]] <- GBS.array@true_parents[[i]]$hap1 + GBS.array@true_parents[[i]]$hap2
    }
    return(GBS.array)
}

out <- data.frame()
for(SIZE in 3:3){
    GBS.array <- sim.array(size.array=SIZE, numloci=10, hom.error = 0.02, het.error = 0.8,
                           rec = 0.25, selfing = 0, imiss = 0.5, misscode = 3)
    GBS.array <- true_other_parents(GBS.array)
    
    inferred_geno_likes <- impute_parent(GBS.array, major.error=0.02, het.error=0.8, minor.error=0.02)
    res <- parentgeno(inferred_geno_likes, oddratio=0.6931472, returnall=TRUE)
    res$true_parent <- GBS.array@true_parents[[SIZE]]$hap1 + GBS.array@true_parents[[SIZE]]$hap2
    
    tem <- data.frame(size=SIZE, error=sum(res$gmax != res$true_parent)) 
    out <- rbind(out, tem)
}


locus=1
outg <- data.frame()
for(i in 1:length(GBS.array@gbs_kids)){
    g <- GBS.array@gbs_parents[[i]][[locus]]
    k <- GBS.array@gbs_kids[[i]][[locus]]
    outg <- rbind(outg, data.frame(otherp=g, kid=k))        
}
i <- i+1
outg$p <- GBS.array@gbs_parents[[i]][[locus]]
GBS.array@snpinfo$frq[[locus]]
