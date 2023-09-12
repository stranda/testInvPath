### plot a few combinations of principle components
### with the observed indicated
#' @title plotPCA
#' @description Plot the results of a PCA.
#' @param pca.ref A reference for the PCA results. TODO: Add more detailed description.
#' @param prop_points Proportion of points to be displayed. Default is 0.2.
#' @param mod_color Color for the model. Default is c(2:5).
#' @return A plot visualizing the PCA results. TODO: Add more detailed description.
#' @export
#' @examples
#' # TODO: Add example usage for plotPCA

plotPCA = function(pca.ref, prop_points=0.2, mod_color=c(2:5), obs_color=6,
                   models=c("Ship30","Ship60","gigasAdmix","gigasML"))
{
    
    obs=pca.ref$obs
    ref=pca.ref$ref
    params=pca.ref$params

    selvec = sample(1:nrow(ref),round(prop_points*nrow(ref)),replace=F)
    params=params[selvec,]
    ref=ref[selvec,]
    
    axes=unique(t(apply(expand.grid(x=1:4,y=1:4),1,sort)))
    axes=axes[axes[,1]!=axes[,2],]

    layout(matrix(c(1, 2), ncol=2), widths=c(18, 4))
    for (i in 1:nrow(axes))
    {
        par(mai=c(0.9,0.9,0.5,0.1))
        
        plot(x=ref[,axes[i,1]],y=ref[,axes[i,2]],
             type="p",pch=c(16,17)[params$SimulIntro+1],col=mod_color[params$introModel],
             main=paste0("PC-",axes[i,1]," on X axis, PC-",axes[i,2]," on Y"),
             xlab=paste0("PC-",axes[i,1]),ylab=paste0("PC-",axes[i,2]),
             xlim=range(c(obs[axes[i,1]],ref[,axes[i,1]])),
             ylim=range(c(obs[axes[i,2]],ref[,axes[i,2]])))
        points(x=obs[axes[i,1]],y=obs[axes[i,2]],cex=2,pch=17,col=obs_color)
        par(mai=c(0.5,0.1,0.5,0.1))
        plot(x=rep(1,length(unique(params$introModel))),
             y=1:4, pch=c(19), cex=16, col=mod_color[1:4],axes=F,xlab="",ylab="",
             main="IntroModel")
        text(x=rep(1,length(unique(params$introModel))),
             y=1:4,models,cex=1.1)
    }
    
}
