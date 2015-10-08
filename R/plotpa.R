plotselfer <- function(sim, kids=1:10, snps=2:99, cols=c("red", "blue"), ...){
    ### input from sim objects returned from SimSelfer
    simp <- sim[[1]]
    simp$id <- 1:nrow(simp)
    simp <- subset(simp, id %in% snps)
    # It sets up the plot but doesn't actually plotting anything
    plot(x=range(simp$id), y=c(1, 5*(length(kids)+1)+2), type = "n", 
         xlab="", ylab="", yaxt="n", ...)
    
    addsnps_mom(genotab=simp[-nrow(simp),], yl=5*(length(kids)+1)+1,  cols=cols)
    #abline(h=4*kids+1, lwd=2, col="grey")
    
    simk <- sim[[2]][kids]
    
    for(i in 1:length(simk)){
        kidi <- simk[[i]]
        infoi <- kidi[[1]]
        
        ### subsetting the selected part of information
        for(t in 1:2){
            for(k in 1:nrow(infoi[[t]])){
                if(infoi[[t]][k, ]$start <= snps[1] & (infoi[[t]][k, ]$end -1) >= snps[1]){
                    idx1 <- k
                    infoi[[t]][k, ]$start <- snps[1]
                }
                if(infoi[[t]][k, ]$start <= snps[length(snps)] & (infoi[[t]][k, ]$end-1) >= snps[length(snps)]){
                    idx2 <- k
                    infoi[[t]][k, ]$end <- snps[length(snps)]
                } 
            }
            infoi[[t]] <- infoi[[t]][idx1:idx2,]
        }
        
        genoi <- kidi[[2]]
        genoi$id <- 1:nrow(genoi)
        genoi <- subset(genoi, id %in% snps)
        
        addsnps_kids(genotab=genoi, yl=5*(length(kids)+1)-5*i, infoi=infoi,  cols=cols)
        
    }   
}


addsnps_mom <- function(genotab=simp, yl=11,  cols=c("red", "blue")){
    # This adds labels to positions (x,y)
    text(x=genotab$id, y=yl, labels=genotab$hap1, col= cols[1]) 
    text(x=genotab$id, y=yl-1, labels=genotab$hap2, col= cols[2]) 
    text(x=genotab$id, y= yl-2, labels=genotab$obs, col="black") 
    sub <- subset(genotab, obs != (hap1+hap2))
    if(nrow(sub) > 0){
        text(x=sub$id, y= yl-2, labels=sub$obs, col="red") 
    }
    abline(h=yl-4, lwd=2, col="grey")
}

addsnps_kids <- function(genotab=simk, yl=11, infoi=infoi,  cols=c("red", "blue")){
    
    hap1 <- infoi[[1]]
    for(j in 1:nrow(hap1)){
        text(x=hap1$start[j]:(hap1$end[j]-1), y=yl, 
             labels= genotab[which(genotab$id == hap1$start[j]):(which(genotab$id == hap1$end[j]) -1),1], 
             col= cols[hap1$hap[j]]) 
    }
    
    hap2 <- infoi[[2]]
    for(k in 1:nrow(hap2)){
        text(x=hap2$start[k]:(hap2$end[k]-1), y=yl-1,
             labels= genotab[which(genotab$id == hap2$start[k]):(which(genotab$id == hap2$end[k]) -1),2],
             col= cols[hap2$hap[k]]) 
    }
    
    genotab <- genotab[-nrow(genotab),]
    text(x=genotab$id, y=yl-2, labels=genotab$obs, col= "black") 
    sub <- subset(genotab, hap1+hap2 != obs)
    if(nrow(sub) >0){
        text(x=sub$id, y=yl-2, labels=sub$obs, col= "red") 
    }
    
    abline(h=yl-4, lwd=2, col="grey")
}