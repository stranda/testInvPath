######
#####
####
## This is a script that runs a particular number of reps of coalescent simulation and
## the calculation of summary statistics.
##

args <- commandArgs(trailingOnly = TRUE)
### Check if an argument is provided, and if it's numeric
if (length(args) == 2 && !is.na(as.numeric(args[2]))) {
    cores <- as.numeric(args[1])
    seed <- round(as.numeric(args[2]))
} else {
    cores <- 1
    slurmArrayID = Sys.getenv("SLURM_ARRAY_TASK_ID")
    if (slurmArrayID!="")
    {
        slurmJobID=Sys.getenv("SLURM_JOB_ID")
        seed = as.integer(paste0(slurmArrayID,slurmJobID))
    }
    else seed <- NULL
}
####set the rng seed based on either the command line 2nd param or NULL
print("seed:")
print(seed)

set.seed(seed)


library(parallel,quiet=T)
library(testInvPath,quiet=T)
library(strataG,quiet=T)
library(phyclust)

source("datafiles.R") #defines mname, genofile, sources and intros

#first set up the basic parameters for these simulations


params = setupReps(mname=paste0(abspath,"/",mname),
                   genofile=paste0(abspath,"/",genofile),
                   native_range_topology=nativeTopology, #present in datafiles.R
                   mainparams(cores=cores,
                              nreps=50,
                              demeN=20,
                              demeI=20,
                              sources=sources,
                              intros=intros,
                              dataType=dataType,
                              species=species,
                              fsc_exec=fsc_exec,
                              exclude="nonSource",
                              popPairwise=popPairwise
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
        runout=fscRun(simp,exec=fsc_exec,dna.to.snp=F,seed=round(runif(1,1,10000000)),tree=use.seqgen)
        if (use.seqgen)
        {
            tree=read.nexus(paste0(simp$folder,"/strataG.fsc/strataG.fsc_1_true_trees.trees"))
            res=seqgen2gtype(tree,params,demes=simp$settings$demes,temp.file=paste0(simp$folder,"/tmpseq.phylip"),mu=round(10*priors$mut_rate,8))
        } else {
            res=fsc2gtypes(runout,marker='dna',as.genotypes=F,concat.dna=F)
            res=res[,1,] #only have 1 sequence in emprirical data
        }
        
    }   else if (tolower(params$dataType)=="snp") {
        runout=fscRun(simp,exec=fsc_exec,max.snps=1000,dna.to.snp=T,seed=round(runif(1,1,10000000)),tree=T)
        res=fsc2gtypes(runout,marker="snp")
    } else if (tolower(params$dataType)=="microsatellite") {
        runout=fscRun(simp,exec=fsc_exec,seed=round(runif(1,1,10000000)))
        res=fsc2gtypes(runout,marker="microsat")
    }  else stop("incorrect data type specified in runreps")

    fscCleanup(runout$label,runout$folder) #remove fsc files
    elapsed=Sys.time()-strt
    print(c(elapsed,unlist(priors)))
    c(unlist(priors),summary_stats(res,params$meta,dataType,popPairwise=params$popPairwise))
}))

write.table(file=paste0("reference",round(runif(1,min=0,max=1000000)),".csv"),row.names=F,sep=",",data.frame(rdf))
