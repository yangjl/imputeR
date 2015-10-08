### make a transition matrix
MakeTransMat = function(stayPrb)
{
    nState = length(stayPrb)
    transMat = matrix(nrow=nState, ncol=nState, data=0)
    for(k in 1:nState){
        transMat[k,] = (1-stayPrb[k]) / (nState - 1)
        transMat[k,k] = stayPrb[k]
    }
    return( transMat )
}
### generate sequence
GenSeq = function(transMat, seqLen, seqMeans, seqSDs = NULL, initDistn = NULL)
{
    nState = ncol(transMat)
    if(is.null(initDistn)) initDistn = rep(1/nState, nState)
    if(is.null(seqSDs)) seqSDs = rep(1, nState)
    rQ = rep(0, seqLen)
    rO = rep(0, seqLen)
    rQ[1] = sample(1:nState, 1, prob=initDistn)
    for(j in 2:seqLen){
        rQ[j] = sample(1:nState, 1, prob=transMat[rQ[j-1],] )
    }
    for(j in 1:nState){
        i1 = which(rQ == j)
        if(length(i1) > 0){
            rO[i1] = rnorm( length(i1), mean = seqMeans[j], sd=seqSDs[j] )
        }
    }
    return( list(O=rO, Q=rQ) )
}

#fwdPrbs = ComputeFwd(x$O, transMat = tM, seqMeans = seqMeans, seqSDs = seqSDs)
ComputeFwd = function(obs, transMat, seqMeans, seqSDs = NULL, initDistn = NULL)
{
    nState = ncol(transMat)
    seqLen = length(obs)
    if(is.null(initDistn)) initDistn = rep(1/nState, nState)
    if(is.null(seqSDs)) seqSDs = rep(1, nState)
    rAlpha = matrix(nrow=nState, ncol=seqLen, data=0)
    rEmitPrbs = matrix(nrow=nState, ncol=seqLen, data=0)
    rEmitPrbs[,1] = rAlpha[,1] = dnorm(obs[1], mean=seqMeans, sd=seqSDs)
    for(j in 2:seqLen){
        rEmitPrbs[,j] = dnorm(obs[j], mean=seqMeans, sd=seqSDs)
        rAlpha[,j] = rEmitPrbs[,j] * (transMat %*% rAlpha[,j-1])
    }
    return( list(emitPrbs=rEmitPrbs, fwdProbs=rAlpha) )
}
SampleSeq = function(nSamp, transMat, fwdProbs)
{
    seqLen = ncol(fwdProbs)
    nState = ncol(transMat)
    if(nSamp > 1){
        rQ = matrix(nrow=seqLen, ncol=nSamp, data=0)
        for(j in 1:nSamp){
            rQ[,j] = SampleSeq(1, transMat, fwdProbs)
        }
        return( rQ )
    }else{
        a = fwdProbs
        nState = ncol(transMat)
        rQ = rep(0, seqLen)
        pr1 = a[,seqLen]/sum(a[,seqLen])
        rQ[seqLen] = sample(1:nState, 1, prob=pr1 )
        for(j in rev(1:(seqLen-1)) ){
            pr1 = a[,j] * transMat[ , rQ[j+1] ]
            pr1 = pr1 / sum(pr1)
            rQ[j] = sample(1:nState, 1, prob=pr1 )
        }
        return( rQ )
    }
}

