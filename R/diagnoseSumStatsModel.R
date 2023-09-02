#' @title obsModelDist
#' @description Calculate the distance between observed and model distributions.
#' @param refObj An object containing reference data. TODO: Add more detailed description.
#' @param plot Logical. If TRUE, plots the differences. Default is FALSE.
#' @return A vector of distances. TODO: Add more detailed description.
#' @export
obsModelDist = function(refObj, plot=FALSE)
{
    ref=refObj$ref
    params=refObj$params
    obs=refObj$obs

    if (!refObj$scale.var)
    {
        rscale=apply(rbind(obs,ref),2,scale)
        obs=data.frame(t(rscale[1,]))
        ref=data.frame(rscale[-1,])
        names(obs)=names(ref)
        rm(rscale)
    }
    dists=colMeans(ref)-obs
    if (plot)
    {
        tmp=c(unlist(dists[1,]))
        names(tmp)=names(dists)
        tmp=abs(tmp)
        tmp=sort(tmp,decreasing=T)
        par(mfrow=c(1,2))
        dotchart(tmp[1:20],main="Absolute difference between\n mean sims and obs",
                  xlim=range(tmp[1:40]),xlab="z-score diff")
        dotchart(tmp[21:40],main="Absolute difference between\n mean sims and obs",
                 xlim=range(tmp[1:40]),xlab="z-score diff")
        par(mfrow=c(1,1))
    }
    
    dists
}

