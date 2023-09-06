## build a reference table by pca transforming summary stats
##
#' @title make_ref_table
#' @description Create a reference table based on genetic data. TODO: Add more detailed description.
#' @param genofile Path to the genetic data file. TODO: Add more detailed description.
#' @param mname TODO: Describe this parameter.
#' @param intros TODO: Describe this parameter.
#' @param sources TODO: Describe this parameter.
#' @param csvpath Path to save the CSV output. Default is '.'.
#' @param dataType Type of genetic data. Default is 'microsatellite'.
#' @param na.prop Proportion of NAs. Default is 0.01.
#' @param scale.var Logical. Whether to scale the variables. Default is TRUE.
#' @param verb       Logical. Whether to print details (default FALSE)
#' @return A reference table based on genetic data. TODO: Add more detailed description.
#' @export
#' @examples
#' # TODO: Add example usage for make_ref_table

make_ref_table = function(genofile, mname, intros, sources, csvpath=".",
                          dataType="microsatellite",
                          na.prop=0.01,scale.var=T,
                          verb=F)
{
    ## calculate the observed summary stats
    meta1=read.csv(mname)
    meta1$longpop=meta1$pop
    gin=readRDS(genofile)
    strataG::setOther(gin,"popnames")=strataG::getStrata(gin)
    gin.orig=gin
    rmpops = c(meta1$pop[!meta1$source%in%c(sources,intros)])
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
    meta$intro = (meta$source %in% intros) 
    meta = meta[,c("longpop","intro","source","idnum","pop")]


    obs=summary_stats(gin,meta,dataType)

    

    ## read in simulated summary stats

    files = list.files(path=csvpath,pattern="*.csv")

    dat=NULL
    for (f in files)
    {
        print(f)
        dat=rbind(dat,data.table::fread(paste0(csvpath,'/',f)))
    }

    dat=as.data.frame(dat)
    params=dat[,1:(ncol(dat)-length(obs))]
    ref=dat[,(ncol(params)+1):ncol(dat)]

    ## remove invariant columns from ref
    invCol=which(apply(ref,2,function(x){var(x)})==0)
    ref=ref[,-invCol]
    obs=obs[-invCol]

    
    ##remove columns that have lots of NAs (or any obs NAs) #same for Inf
    naCount = apply(ref,2,function(x){sum((!is.finite(x)))})
    badcols = unique(c(which((naCount/nrow(ref))>na.prop),
                       which(!is.finite(obs))))
    ref = data.frame(ref) [,-1*badcols]
    obs = data.frame(t(obs [-1*badcols]))
    
    ##remove rows that have any NAs (after removing bad columns above)
    ## have to remember to remove same rows from params as well
    cc = complete.cases(ref)
    ninf=colSums(!apply(ref,1,is.finite))==0
    k = cc&ninf
    ref=ref[k,]
    params=params[k,]  

    rm(dat,badcols,cc,naCount)
    if (verb)
    {
        print(summary(ref))
        print(obs)
        print(summary(params)
    }
    
    nobs=obs
    if (scale.var)
    {
        rscale=apply(rbind(obs,ref),2,scale)
        #check that scaling didnt introduce Na
        naCol = apply(rscale,2,function(x){sum((!is.finite(x)))})>0 
        rscale=rscale[,!naCol]
        nobs=data.frame(t(rscale[1,]))
        ref=data.frame(rscale[-1,])
        rm(rscale)
    }

    
    names(nobs)=names(ref)
    list(obs=nobs, params=data.frame(params), ref=ref, scale.var=scale.var)

}

