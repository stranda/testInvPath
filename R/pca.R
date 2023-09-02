## build a reference table by pca transforming summary stats
##
#' @title make_pca_ref_table
#' @description Perform PCA on reference data and return a table.
#' @param refObj An object containing reference data. TODO: Add more detailed description.
#' @param prop.variation Proportion of variation to be captured. Default is 0.99.
#' @param scale.var Logical. Whether to scale the variables. Default is TRUE.
#' @return A table resulting from the PCA. TODO: Add more detailed description.
#' @export
#' @examples
#' # TODO: Add example usage for make_pca_ref_table
make_pca_ref_table = function(refObj, prop.variation=0.99, scale.var=T)
{
    ref=refObj$ref
    params=refObj$params
    obs=refObj$obs
    scale.var = ifelse(refObj$scale.var,F,scale.var) #if already scaled set scale.vat to F otherwise use parameter passed.
    comb=rbind(obs,ref)
    pmod=prcomp(comb,scale=scale.var)
    num.to.keep=which(summary(pmod)$importance[3,]>prop.variation)[1]
    tran=predict(pmod)[,1:num.to.keep]
    obs=data.frame(matrix(tran[1,],nrow=1)) #the usual problems with 1-row rectangular structures
    names(obs)=colnames(tran)
    
    list(obs=obs,ref=as.data.frame(tran[-1,]),params=params,scale.var=scale.var,prop.variation=prop.variation)
}

