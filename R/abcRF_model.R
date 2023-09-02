#' ref2RF Utility Function
#'
#' This internal function converts reference objects for Random Forest.
#'
#' @param refobj A reference object.
#' @param modelColumn A string specifying the model column.
#'
#' @return A data frame prepared for Random Forest.
#' @keywords internal
ref2RF = function(refobj,modelColumn)
{
    ref=data.frame(cbind(refobj$params[,modelColumn],refobj$ref))
    names(ref)[1] = "modelIndex"
    ref$modelIndex=factor(ref$modelIndex)
    ref
}

#' Build an ABC Random Forest Model
#'
#' This function builds an Approximate Bayesian Computation (ABC) Random Forest model.
#'
#' @param refobj A reference object.
#' @param modelColumn A string specifying the model column. Default is "introModel".
#' @param cores An integer specifying the number of cores to use. Default is 1.
#' @param ... Additional arguments.
#'
#' @return An ABC Random Forest model.
#' @export
modelRF = function(refobj,
                   modelColumn="introModel",
                   cores=1,...)
{
    ref = ref2RF(refobj, modelColumn)
    abcrf(modelIndex~.,data=ref,
          group=list(shipping=c("1","2"),gigas=c("3","4")),
          paral=ifelse(cores>1,T,F),ncores=cores)
}

#' Predict Using ABC Random Forest Model
#'
#' This function makes predictions using an ABC Random Forest model.
#'
#' @param rfFit A fitted ABC Random Forest model.
#' @param refobj A reference object.
#' @param modelColumn A string specifying the model column. Default is "introModel".
#' @param cores An integer specifying the number of cores to use. Default is 1.
#'
#' @return Predictions from the ABC Random Forest model.
#' @export
predRFmodel = function(rfFit, refobj, modelColumn="introModel",cores=1)
{
    ref = ref2RF(refobj, modelColumn)
    predict(rfFit,obs=refobj$obs,ref,
            paral=ifelse(cores>1,T,F),ncores=cores,
            paral.predict=ifelse(cores>1,T,F),ncores.predict=cores,
            ntree=5000)
}

#' Plot ABC Random Forest Model
#'
#' This function plots an ABC Random Forest model.
#'
#' @param rfFit A fitted ABC Random Forest model.
#' @param refobj A reference object.
#' @param modelColumn A string specifying the model column. Default is "introModel".
#'
#' @return Plots related to the ABC Random Forest model.
#' @export
plotRFmodel = function(rfFit, refobj, modelColumn="introModel")
{
    ref = ref2RF(refobj, modelColumn)
    plot(rfFit,ref,obs=refobj$obs,pdf=F,nvar=40)
}

#' Calculate Error for ABC Random Forest Model
#'
#' This function calculates the error of an ABC Random Forest model.
#'
#' @param rfFit A fitted ABC Random Forest model.
#' @param refobj A reference object.
#' @param modelColumn A string specifying the model column. Default is "introModel".
#' @param cores An integer specifying the number of cores to use. Default is 1.
#'
#' @return Error metrics for the ABC Random Forest model.
#' @export
errRFmodel = function(rfFit, refobj, modelColumn="introModel",cores=1)
{
    ref = ref2RF(refobj, modelColumn)
    err.abcrf(rfFit,ref, paral=ifelse(cores>1,T,F),ncores=cores )
}
