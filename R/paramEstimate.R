#' paramPosteriors function
#'
#' This function computes the posterior distribution of parameters using Approximate Bayesian Computation (ABC).
#'
#' @param SimObj A simulation object containing reference data, parameters, and observed data.
#' @param cores Number of cores to use for parallel computation. Default is 1.
#' @param method Method to use for ABC. Default is "neuralnet".
#' @param tol Tolerance levels for ABC. Default is c(0.1,0.05,0.01,0.005,0.001,0.005).
#'
#' @return A list of ABC results for each combination of method and tolerance level.
#' Each element of the list is named by concatenating the method and tolerance level with an underscore.
#'
#' @examples
#' # Assuming appropriate data and functions are available:
#' result <- paramPosteriors(SimObj)
#'
#' @export
paramPosteriors = function(SimObj,cores=1,method="neuralnet",
                           tol=c(0.1,0.05,0.01,0.005,0.001,0.005))
{
    ref=SimObj$ref
    params=SimObj$params
    goodCol=apply(params,2,function(x)var(x))>0
    params=params[,goodCol]
    obs=SimObj$obs
    d=expand.grid(method=method,tol=tol)
    
    ret = lapply(1:nrow(d),function(i)
    {
        abc(target=obs,param=params,sumstat=ref,method=method[i],tol=tol[i])
    })
    names(ret)=paste(d$method,d$tol,sep="_")
    ret    
}

