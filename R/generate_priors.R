#' defines prior distributions and generates a pull for each call
#' @description The parameters that are estimated are each chosen from prior distributions in this function.
#' @return a list of prior values
#' @export
getPriors <- function()
{
    tIntro  = runif(1,2,200)
    t0      = runif(1,tIntro+1,10000)
    t1      = runif(1,t0+10,50000)
    t2      = runif(1,t1+1,50000)
    t3      = runif(1,max(c(t1,t2))+1,100000)

    list(ancestralNe = runif(1,1,400), #this is actually a coeffient to multiply nativeNe in the basal event
         introNe = runif(1,5,5000), #effective size of introduced pops
         nativeNe = runif(1,1000,7000), #effective size of native pops
         tIntro=tIntro,t0=t0,t1=t1,t2=t2,t3=t3,
         IntroGrow=0, #-1 * runif(1,0,0.1),
         mut_rate = runif(1,0,0.005),
         introM = runif(1,0,0.01),
         nativeM = runif(1,0,0.01),
         introModel = sample(1:4,1),
         SimulIntro =  sample(c(TRUE,FALSE),1)
        )
}
