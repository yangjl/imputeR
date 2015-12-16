# simulaltion 

true_other_parents <- function(GBS.array){
    for(i in 1:(length(GBS.array@gbs_parents)-1)){
        GBS.array@gbs_parents[[i]] <- GBS.array@true_parents[[i]]$hap1 + GBS.array@true_parents[[i]]$hap2
    }
    return(GBS.array)
}

sim_ip <- function(numloci=1000, selfrate=1, truep=1, outfile="cache/out1.csv"){

    out <- data.frame()
    for(SIZE in 1:100){
        GBS.array <- sim.array(size.array=SIZE, numloci, hom.error = 0.02, het.error = 0.8, selfing=selfrate,
                               rec = 0.25, imiss = 0.5, misscode = 3)
        GBS.array <- true_other_parents(GBS.array)
        GBS.array@pedigree$true_p <- truep
        
        inferred_geno_likes <- impute_parent(GBS.array, major.error=0.02, het.error=0.8, minor.error=0.02)
        res <- parentgeno(inferred_geno_likes, oddratio=0.6931472, returnall=TRUE)
        res$true_parent <- GBS.array@true_parents[[SIZE+1]]$hap1 + GBS.array@true_parents[[SIZE+1]]$hap2
        
        tem <- data.frame(size=SIZE, error=sum(res$gmax != res$true_parent)) 
        out <- rbind(out, tem)
    }
    write.table(out, outfile, sep=",", row.names=FALSE, quote=FALSE)
    return(out)
}

### completely selfed kids
out1 <- sim_ip(numloci=100, selfrate=1, truep=0, outfile="test/simip_out1.csv")

### half selfed kids and half oc (parents unknow)
out2 <- sim_ip(numloci=100, selfrate=0.5, truep=0, outfile="test/simip_out2.csv")

### half selfed kids and half oc (know parents)
out3 <- sim_ip(numloci=100, selfrate=0.5, truep=1, outfile="test/simip_out3.csv")

### complete oc (unknow parents)
out4 <- sim_ip(numloci=100, selfrate=0, truep=0, outfile="test/simip_out4.csv")

### complete oc (know parents)
out5 <- sim_ip(numloci=100, selfrate=0, truep=1, outfile="test/simip_out5.csv")


par(mfrow=c(2,3))
plot(out1[,2])
plot(out2[,2])
plot(out3[,2])
plot(out4[,2])
plot(out5[,2])


