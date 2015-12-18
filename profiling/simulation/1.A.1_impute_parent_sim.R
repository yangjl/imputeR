# simulaltion 

sim_ip <- function(numloci=1000, selfrate=1, outfile=NULL){
    
    
    perr <- gen_error_mat(major.error=0.02, het.error=0.8, minor.error=0.02)
    kerr <- perr
    
    out <- data.frame()
    for(SIZE in 1:100){
        GBS.array <- sim.array(size.array=SIZE, numloci, hom.error = 0.02, het.error = 0.8, selfing=selfrate,
                               rec = 0.25, imiss = 0.5, misscode = 3)
        #GBS.array <- true_other_parents(GBS.array)
        #GBS.array@pedigree$true_p <- truep
        
        inferred_geno_likes <- impute_parent(GBS.array, perr, kerr)
        res <- parentgeno(inferred_geno_likes, oddratio=0.6931472, returnall=TRUE)
        res$true_parent <- GBS.array@true_parents[[SIZE+1]]$hap1 + GBS.array@true_parents[[SIZE+1]]$hap2
        
        tem <- data.frame(size=SIZE, error=sum(res$gmax != res$true_parent)) 
        out <- rbind(out, tem)
    }
    if(!is.null(outfile)){
        write.table(out, outfile, sep=",", row.names=FALSE, quote=FALSE)
    }
    
    return(out)
}

### completely selfed kids
out1 <- sim_ip(numloci=100, selfrate=1, outfile="test/simip_out1.csv")

### half selfed kids and half oc (parents unknow)
out2 <- sim_ip(numloci=100, selfrate=0.5, outfile="test/simip_out2.csv")

### complete oc (unknow parents)
out3 <- sim_ip(numloci=100, selfrate=0, outfile=NULL)


par(mfrow=c(1,3))
plot(out1[,2]/100, main="Parental Imputation", xlab="Family Size", ylab="Imputation Error Rate")
plot(out2[,2]/100, main="Parental Imputation", xlab="Family Size", ylab="Imputation Error Rate")
plot(out3[,2]/100, main="Parental Imputation", xlab="Family Size", ylab="Imputation Error Rate")

