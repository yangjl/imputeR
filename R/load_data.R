#######################################
maf_missing <- function(wgs, gbs){
    lmiss1 <- apply(wgs[, 9:27], 1, function(x) return(sum(x==3)/19))
    lmiss2 <- apply(gbs[, 3:21], 1, function(x) return(sum(x==3)/19))
    
    imiss1 <- apply(wgs[, 9:27], 2, function(x) return(sum(x==3)/301249))
    imiss2 <- apply(gbs[, 3:21], 2, function(x) return(sum(x==3)/301249))
    
    getmaf <- function(dmx){
        unlist(apply(dmx, 1, function(x){
            x <- as.numeric(as.character(x))
            x <- x[x!=3]
            if(length(x) >0 ){
                c0 <- sum(x == 0)
                c1 <- sum(x == 1)
                c2 <- sum(x == 2)
                return(min(c(2*c0+c1, c1+2*c2))/(2*(c0 + c1 + c2)) )
            } }))
    }
    
    maf1 <- getmaf(wgs[, 9:27])
    maf2 <- getmaf(gbs[, 3:21])
    
    outfile="cache/teo_gbs_wgs.RData"
    message(sprintf("###>>> Data write to: [ %s]", outfile))
    save(file=outfile, list=c("lmiss1", "lmiss2", "imiss1", "imiss2", "maf1", "maf2"))    
    
}


###########################################################
comp_alleles <- function(wgs, gbs){
    
    wgs <- wgs[order(wgs$snpid2), ]
    gbs <- gbs[order(gbs$snpid2), ]
    gbs$snpid2 <- gsub("S", "", gbs$snpid2)
    
    nms <- names(gbs)[-1:-2]
    
    heterr <- hettot <- homerr <- homtot <- 0
    for(i in 1:length(nms)){
        out <- merge(wgs[, c("snpid2", nms[i]) ], gbs[, c("snpid2", nms[i]) ], by="snpid2")
        names(out) <- c("snpid", "g1", "g2")
        
        out <- subset(out, g1 !=3 & g2 != 3)
        
        if(nrow(out) >0){
            heterr <- heterr + nrow(subset(out, g1 == 1 & g1 != g2))
            hettot <- hettot + nrow(subset(out, g1 == 1))
            homerr <- homerr + nrow(subset(out, g1 !=1 & g1 != g2))     
            homtot <- homtot + nrow(subset(out, g1 !=1))
        }
        
    }
    message(sprintf("###>>> Heterozygote error rate [ %s ] and Homozygote error rate [ %s ]",
                    round(heterr/hettot, 3)*100, round(homerr/homtot, 3)*100))
    message(sprintf("###>>> het err=[ %s ]; het tot=[ %s ]; hom err=[ %s ]; hom err=[ %s ]", 
                    heterr, hettot, homerr, homtot))
    return(c(heterr, hettot, homerr, homtot))
    
}

##################################################################
recode <- function(){
    ob <- load("largedata/wgs_teo19.RData")
    ### steo: 396818; v
    info <- read.csv("largedata//teo_info.csv")
    info$snpid <- gsub("S", "", info$snpid)
    info <- merge(info, v[, 5:6], by.x="snpid", by.y="snpid2")
    names(info)[1] <- "snpid2"
    comp <- merge(steo[, c("snpid", "major", "minor")], info[, c(11, 1:3)], by.x="snpid", by.y="snpid3")
    
    message(sprintf("###>>> WGS [ %s ] | GBS [ %s ] | shared [ %s ]", nrow(steo), nrow(info), nrow(comp)))
    
    ### Teo19 WGS V3 and V4 are major/minor
    idx <- which((comp$major == comp$ref & comp$minor == comp$alt) | (comp$major == comp$alt & comp$minor == comp$ref))
    message(sprintf("###>>> consistent SNP calling [ %s ]", length(idx)))
    
    steo <- merge(comp[idx, c(1,4:6)], steo, by="snpid")
    ### recoding ATCG=> 0, 1, 2
    for(i in 9:ncol(steo)){
        steo[, i] <- as.character(steo[, i])
        steo$a1 <- gsub(".$", "", steo[, i])
        steo$a2 <- gsub("^.", "", steo[, i])
        
        steo[steo[, i]!= "NN" & steo$a1 == steo$alt, ]$a1 <- 1
        steo[steo[, i]!= "NN" & steo$a1 == steo$ref, ]$a1 <- 0
        steo[steo[, i]!= "NN" & steo$a2 == steo$alt, ]$a2 <- 1
        steo[steo[, i]!= "NN" & steo$a2 == steo$ref, ]$a2 <- 0
        
        steo[steo$a1 == "N", ]$a1 <- 1.5
        steo[steo$a2 == "N", ]$a2 <- 1.5
        
        steo[, i] <- as.numeric(as.character(steo$a1)) + as.numeric(as.character(steo$a2))
    }
    steo$snpid <- paste0("S", steo$snpid)
    return(steo)
}

################
gbsgeno <- function(steo){
    
    ### SNP matrix comparison
    library(parallel)
    library(devtools)
    options(mc.cores=NULL)
    load_all("~/bin/tasselr")
    load_all("~/bin/ProgenyArray")
    
    ob2 <- load("largedata/cj_data.Rdata")
    genos <- geno(teo)
    
    nms <- gsub("_1\\:.*|_mrg\\:.*", "", colnames(genos))
    subgeno <- genos[, which(nms %in% names(steo)[9:27])]
    subgeno[is.na(subgeno)] <- 3
    subgeno <- as.data.frame(subgeno)
    names(subgeno) <- gsub("_1\\:.*|_mrg\\:.*", "", names(subgeno))
    message(sprintf("###>>> GBS of [ %s ] SNPs and [ %s ] plants", nrow(subgeno), ncol(subgeno)))
    
    subgeno$snpid2 <- as.character(row.names(subgeno))
    steo$snpid2 <- paste0("S", steo$snpid2)
    tem <- merge(steo[, 1:2], subgeno, by = "snpid2")
    message(sprintf("###>>> Common SNPs [ %s ] ", nrow(tem)))
    return(tem)
}


