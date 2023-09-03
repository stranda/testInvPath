#' Calculate Posterior Probabilities for Different Models
#'
#' This function computes the posterior probabilities of different models 
#' given observed data using Approximate Bayesian Computation (ABC).
#'
#' @param refobj A reference object containing the observed data and reference table.
#' @param method A character vector indicating the regression method(s). Default is "neuralnet".
#' @param tol A numeric vector indicating the tolerance levels for ABC. 
#'            Default is c(0.1,0.05,0.01,0.005,0.001,0.005).
#' @param modelColumn A string specifying the model column. Default is "introModel".
#' @param cores An integer specifying the number of cores to use for parallel computation. Default is 1.
#' @param ... Additional arguments passed to the regression methods.
#'
#' @return A list containing the posterior probabilities for different models under 
#'         various method and tolerance settings.
#' @export

ModelPostProbs = function(refobj,
                          method=c("neuralnet"),
                          tol=c(0.1,0.05,0.01,0.005,0.001,0.005),
                          modelColumn="introModel",
                          cores=1,...)
{
    
    d=expand.grid(method=method,tol=tol)
    modvec=refobj$params[,modelColumn]

    modsel = parallel::mclapply(1:nrow(d),mc.cores=cores,function(i)
    {
        if (d$method[i]=="neuralnet")
        {
            if (exists("sizenet")) sn = sizenet else sn = 5
            if (exists("MaxNWts")) mw = MaxNWts else mw = 1000
            if (exists("maxit")) mi = maxit else mi=2500
            
            col=ncol(refobj$ref)
            tv=round(col*sn*1.01)
            if (tv>mw)
            {
                print(paste("There are so many columns the maximume weights are exceeded",tv,"versus",mw))
                print("Adjusting maxweights up, but this will slow analysis")
                mw=tv
            }
            
            postpr(refobj$obs,modvec,refobj$ref,tol=d$tol[i],
                   method=d$method[i],maxit=mi,MaxNWts=mw)
        }
        else
            postpr(refobj$obs,modvec,refobj$ref,tol=d$tol[i],
                   method=d$method[i])

    })
    names(modsel)=paste(d$method,d$tol,sep="_")
    modsel 
}


#' Cross-Validation for Model Choice Using ABC
#'
#' This function performs cross-validation for model choice using the `cv4postpr` function 
#' from the `abc` package to estimate misclassification probabilities for each model.
#'
#' @param refobj A reference object containing the observed data and reference table.
#' @param method A character vector indicating the regression method(s). Default is "neuralnet".
#' @param tol A numeric vector indicating the tolerance levels for ABC. 
#'            Default is c(0.1,0.05,0.01,0.005,0.001,0.005).
#' @param modelColumn A string specifying the model column. Default is "introModel".
#' @param cores An integer specifying the number of cores to use for parallel computation. Default is 1.
#' @param ... Additional arguments passed to the regression methods.
#'
#' @return A list containing cross-validation results for different models under 
#'         various method and tolerance settings.
#' @export

CVModelPostProbs = function(refobj,
                          method=c("neuralnet"),
                          tol=c(0.1,0.05,0.01,0.005,0.001,0.0005),
                          modelColumn="introModel",
                          nval=10,
                          cores=1,...)
{
    
    d=expand.grid(method=method,tol=tol)
    modvec=refobj$params[,modelColumn]

    cv_results = parallel::mclapply(1:nrow(d),mc.cores=cores,function(i)
    {
        if (d$method[i]=="neuralnet")
        {
            if (exists("sizenet")) sn = sizenet else sn = 5
            if (exists("MaxNWts")) mw = MaxNWts else mw = 1000
            if (exists("maxit")) mi = maxit else mi=2500
            
            col=ncol(refobj$ref)
            tv=round(col*sn*1.01)
            if (tv>mw)
            {
                print(paste("There are so many columns the maximum weights are exceeded",tv,"versus",mw))
                print("Adjusting maxweights up, but this will slow analysis")
                mw=tv
            }
            
            cv4postpr(modvec,refobj$ref,tol=d$tol[i],
                   method=d$method[i],nval=nval,maxit=mi,MaxNWts=mw)
        }
        else
            cv4postpr(modvec,refobj$ref,tol=d$tol[i],
                   method=d$method[i],nval=nval)

    })
    names(cv_results)=paste(d$method,d$tol,sep="_")
    cv_results 
}


