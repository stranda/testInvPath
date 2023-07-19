#' setup a single simulation rep for abc
#' 
#' @param mname name of metadata file for pops
#' @param genofile name of gtypes file
#' @param params a list of parameters that come from mainparams()
#'
#' @description This function takes the meta data and the geentic data that are placed
#' in the 'data' subfolder to create a set of parameters that are common to all simulation
#' replicates.  This function only need to be called once when a set of simulations starts. 
#' 
#' @return a list of parameters needed for coalescent simulation creation
#' 
#' @export

setupReps = function(mname,genofile,
                  params = mainparams())
{
    meta1=read.csv(mname)
    meta1$longpop=meta1$pop
    gin=readRDS(genofile)
    strataG::setOther(gin,"popnames")=strataG::getStrata(gin)
    gin.orig=gin
    
    rmpops = c(meta1$pop[!meta1$source%in%c(params$sources,params$intros)])
    
    keepvec = !(strataG::getStrata(gin)%in% rmpops)
    gin=gin[keepvec,,,drop=T]#remove non-source
    meta1=meta1[!(meta1$pop %in% rmpops),]

    
    meta = meta1
    meta$idnum = 0:(nrow(meta)-1)


    pop=strataG::getOther(gin,"popnames")
    pop=pop[pop %in% meta$pop]
    demedf = data.frame(longpop=unique(pop))
    
    samples = as.data.frame(table(pop))
    names(samples) = c("longpop","sample.size")
    meta$intro = (meta$source %in% params$intros) 
    meta = meta[,c("longpop","intro","source","idnum","pop")]

    native_range_topology=nativeHistory(sources=sources)

    introTbls = list(shipTbl30=shipTbl30(intros,sources),
                     shipTbl60=shipTbl60(intros,sources),
                     gigasTblAdmix=gigasTblAdmix(intros,sources),
                     gigasTblML=gigasTblML(intros,sources),
                     equalSourceTbl=equalSourceTbl(intros,sources)
                     )

    simparams = list(gin=gin, species=params$species, meta=meta,introTbls=introTbls,
                     native_range_topology=native_range_topology,demeN=params$demeN,
                     demeI=params$demeI,samples=samples,dataType=params$dataType,
                     fsc_exec=params$fsc_exec,cores=params$cores, nreps=params$nreps)

    
}
