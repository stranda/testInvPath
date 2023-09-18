plotCoalTree = function(tree,demes,meta)
{
    src=meta$source
    names(src)=meta$pop
    tip.color=unclass(factor(src[demes$deme.name[sapply(strsplit(tree$tip.label,"_"),function(x){as.numeric(x[1])})]]))
    plot.phylo(tree,tip.color=tip.color)

    usr <- par("usr")
    x_legend <- usr[1] + 0.10 * (usr[2] - usr[1])
    y_legend <- usr[3] + 0.90 * (usr[4] - usr[3])

    legend(x=x_legend,y=y_legend,
           legend=unique(src[demes$deme.name[sapply(
                                       strsplit(tree$tip.label,"_"),function(x){as.numeric(x[1])})]]),
           lty=1,lwd=7,cex=1.8,
           col=unique(tip.color))

}
