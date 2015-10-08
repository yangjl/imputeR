plotphasing <- function(sim, kids=1:10, snps=2:99, cols=c("red", "blue"), plotphasing=FALSE, phase, ...){
    ### input from sim objects returned from SimSelfer
    simp <- sim[[1]]
    simp$id <- 1:nrow(simp)
    simp <- subset(simp, id %in% snps)
    # It sets up the plot but doesn't actually plotting anything
    plot(x=range(simp$id), y=c(1, 5*(length(kids)+1)+2), type = "n", 
         xlab="", ylab="", yaxt="n", ...)
    
    phase$id <- 1:nrow(phase)
    phase <- subset(phase, id %in% snps)
    phase <- phase[-nrow(phase), ]
    addsnps_mom2(genotab=simp[-nrow(simp),], yl=5*(length(kids)+1)+2,  cols=cols, plotphasing, phase)
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


addsnps_mom2 <- function(genotab=simp, yl=11,  cols=c("red", "blue"), plotphasing, phase){
    # This adds labels to positions (x,y)
    text(x=genotab$id, y=yl+1, labels=genotab$hap1, col= cols[1]) 
    text(x=genotab$id, y=yl, labels=genotab$hap2, col= cols[2]) 
    text(x=genotab$id, y= yl-1, labels=genotab$obs, col="black") 
    sub <- subset(genotab, obs != (hap1+hap2))
    if(nrow(sub) > 0){
        text(x=sub$id, y= yl-1, labels=sub$obs, col="red") 
    }
    if(plotphasing){
        idx <- which.max(c(cor(phase$h1, genotab$hap1), cor(phase$h1, genotab$hap2)))
        text(x=phase$id, y=yl-2, labels=phase$h1, col= cols[idx]) 
        text(x=phase$id, y=yl-3, labels=phase$h2, col= cols[3-idx])
        sub <- subset(phase, h1 != genotab[, idx])
        if(nrow(sub) > 0){
            text(x=sub$id, y=yl-2, labels=sub$h1, col= "yellow")
            text(x=sub$id, y=yl-2, labels=sub$h2, col= "yellow")
        }
    }
    
    abline(h=yl-5, lwd=2, col="grey")
}

