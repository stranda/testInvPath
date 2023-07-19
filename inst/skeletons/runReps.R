######
#####
####
## This is a script that runs a particular number of reps of coalescent simulation and
## the calculation of summary statistics.
##

library(parallel)
library(testInvPath)
library(strataG)

source("datafiles.R") #defines mname, genofile, sources and intros

#first set up the basic parameters for these simulations

params = setupReps(mname=paste0(abspath,"/",mname),
                   genofile=paste0(abspath,"/",genofile),
                   mainparams(cores=1,
                              nreps=1,
                              demeN=0,
                              demeI=0,
                              sources=sources,
                              intros=intros,
                              dataType=dataType,
                              species=species,
                              fsc_exec=fsc_exec,
                              exclude="nonSource"
                              )
                   )

rdf = do.call(rbind,mclapply(1:params$nreps,mc.cores=params$cores,function(i)
{
    priors=getPriors()
    print(priors$introModel)
    strt=Sys.time()
    simp=createCoalObject(params,priors,params$demeN,params$demeI)

    if (tolower(params$dataType)=="sequence")
    {
        runout=fscRun(simp,exec=fsc_exec,dna.to.snp=F)
        res=fsc2gtypes(runout,marker='dna',as.genotypes=F,concat.dna=F)
        res=res[,1,] #only have 1 sequence in emprirical data
    }   else if (tolower(params$dataType)=="snp") {
        runout=fscRun(simp,exec=fsc_exec,max.snps=1000,dna.to.snp=T)
        res=fsc2gtypes(runout,marker="snp") 
    } else if (tolower(params$dataType)=="microsatellite") {
        res=fsc2gtypes(fscRun(simp,exec=fsc_exec),marker="microsat")
    }  else stop("incorrect data type specified in runreps")

    fscCleanup(runout$label,runout$folder) #remove fsc files
    elapsed=Sys.time()-strt
    print(c(elapsed,unlist(priors)))
    c(unlist(priors),summary_stats(res,params$meta,dataType))
}))

write.table(file=paste0("reference",round(runif(1,min=0,max=1000000)),".csv"),row.names=F,sep=",",data.frame(rdf))
