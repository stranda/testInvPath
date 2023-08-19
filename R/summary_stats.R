###
### summary stat calculator  There are separate functions for different dataTypes
###


#' calculate summary stats
#' @param gin is a gtypes object with genetic data
#' @param meta is a dataframe with metadata
#' @param dataType is a character that defines the type of genetic
#'     data
#' @description this function calls different summary stats functions
#'     based on the value of dataType (possibilities include
#'     "sequence", "snp", "microsatellite")
#' @return vector of summary statistics
#' @export
summary_stats=function(gin,meta,dataType)
{
    if (tolower(dataType)=="sequence")
        summary_stats_seq(gin,meta)
    else if (tolower(dataType)=="snp")
        summary_stats_snp(gin,meta)
    else if (tolower(dataType)=="microsatellite")
        summary_stats_microsatellite(gin,meta)
    else stop("incorrect data type specified in summary_stats()")
}

#' calculate summary stats for sequence data (like mtDNA or cpDNA)
#' @param gin is a gtypes object with sequence data
#' @param meta is a dataframe with metadata
#' @description calculates a lot of population genetic summary statistics
summary_stats_seq=function(gin,meta)
{
    ###nucleotide diversity among and within pops
    overallDiv=c(mean(nucleotideDiversity(gin)))
    names(overallDiv)="overallDiversity"
    withinPopDiv=sapply(unique(getStrata(gin)),function(s)
    {
        res=mean(nucleotideDiversity(gin[,,s],simplify=T))
        res
    })
    names(withinPopDiv)=paste0(names(withinPopDiv),".within")
    withinPopDiv[is.na(withinPopDiv)]=0 #replace NA with no diversity
    nd=nucleotideDivergence(gin)
    withinPi = nd$within$mean
    withinPi[is.na(withinPi)]=0 #replace NA with no diversity
    names(withinPi)=nd$within$stratum
    
    amongPi=nd$between$mean
    names(amongPi)=paste0(nd$between$strata.1,".",nd$between$strata.2,".pi")
    
    ##haplotype diversity
    overallHapDiv=pegas::hap.div(getSequences(gin))
    names(overallHapDiv)="overallHapDiv"
    withinPopHapDiv=sapply(unique(getStrata(gin)),function(s)
    {
        ret = 0
        if (!is.na(withinPopDiv[paste0(s,".within")]))
            if (withinPopDiv[paste0(s,".within")]!=0)
                if (table(getStrata(gin))[s]>1)
                    ret =  pegas::hap.div(getSequences(labelHaplotypes(gin[,,s])))
        ret
    })
    names(withinPopHapDiv)=paste0(names(withinPopDiv),".hap")
    overallStruct=overallTest(nosingles(gin),nrep=0)$result[3,1]
    pwStruct=sapply(pairwiseTest(nosingles(gin),nrep=0,quiet=T),function(x){y=x$result[3,1];names(y)=paste(names(x$strata.freq),collapse=".");y})
    names(overallStruct) = "overallPhiST"
    names(pwStruct)=paste0(names(pwStruct),".PhiST")

###source and intro regions
    pops=data.frame(cbind(pop=strataG::getStrata(gin),ind=strataG::getIndNames(gin)))
    smap = meta$source
    names(smap)=meta$pop
    pops$region = smap[pops$pop]
    regions = pops$region
    names(regions)=pops$ind
    strataG::setStrata(gin) = regions #should convert strata to sources
    withinRegionDiv=sapply(unique(getStrata(gin)),function(s)
    {
        res=mean(nucleotideDiversity(gin[,,s],simplify=T))
        res
    })
    names(withinRegionDiv)=paste0(names(withinRegionDiv),".within")
    withinRegionHapDiv=sapply(unique(getStrata(gin)),function(s)
    {
        ret = 0
        if (!is.na(withinRegionDiv[paste0(s,".within")]))
            if (withinRegionDiv[paste0(s,".within")]!=0)
                if (table(getStrata(gin))[s]>1)
                    ret =  pegas::hap.div(getSequences(labelHaplotypes(gin[,,s])))
        ret
    })
    names(withinRegionHapDiv)=paste0(names(withinRegionDiv),".hap")
    pwStruct.Region=sapply(pairwiseTest(nosingles(gin),nrep=0,quiet=T),function(x){y=x$result[3,1];names(y)=paste(names(x$strata.freq),collapse=".");y})
    names(pwStruct.Region)=paste0(names(pwStruct.Region),".PhiST")
        
    r = c( overallDiv,  withinPopDiv, withinPi,  amongPi, overallHapDiv,
          withinPopHapDiv, overallStruct, pwStruct,withinRegionHapDiv, pwStruct.Region)
    r[order(names(r))]
}


###takes a strataG 'gtypes' object and returns one with all strata with a single
###individual removed
nosingles = function(gin)
{
    st = table(getStrata(gin))
    if (min(st)<2)
        gin[,,names(st)[which(st>1)]]
    else
        gin
}


###
### summary stat calculator for snps.  this is a sfs version (less tested)
###

summary_stats_sfs=function(gin,meta)

{
    snps = gin@data %>% group_by(id,locus) %>%
        summarize(snp=sum(as.numeric(allele))) %>%
        pivot_wider(id_cols=1,names_from=2,values_from = 3)
    snps=as.data.frame(snps)
    lp=getStrata(gin)[snps$id]
    snps=data.frame(cbind(id=snps[,1],longpop=unname(unlist(lp)),
                          apply(snps[,-1],2,function(x) if(mean(x)>1) abs(x-2) else x)))
    for (i in 3:ncol(snps)) snps[,i] = as.integer(snps[,i])
    if (sum(apply(snps[,-1:-2],2,function(x) length(unique(x))>1)) != (ncol(snps)-2))
        print("less than the expected polymorphic snps")
    sfs.longpop = sfs(as.data.frame(snps))

    source=unique(meta[,c("source","longpop")])
    svec=source$source
    names(svec)=source$longpop
    tmpvec=snps[,2]
    snps[,2]=svec[snps[,2]]
    names(snps)[2]="src"
    sfs.src=sfs(as.data.frame(snps))
    snps[,2]=tmpvec
    
    intro=unique(meta[,c("intro","longpop")])
    svec=ifelse(intro$intro,"Intro","Native")
    names(svec)=intro$longpop
    tmpvec=snps[,2]
    snps[,2]=svec[snps[,2]]
    names(snps)[2]="status"
    sfs.status=sfs(as.data.frame(snps))
    snps[,2]=tmpvec

###    c(unlist(sfs.longpop),unlist(sfs.src),unlist(sfs.status))
    c(unlist(sfs.src))
}

###
### summary stat calculator for snps
###
summary_stats_trad = function(gin,meta)
{
    meta.orig=meta
    rownames(meta)=meta$longpop
    overallHet = mean(heterozygosity(gin)[,2])
    names(overallHet)="overallHet"
    ostruct=overallTest(gin,nrep=0)$result[,1]
    
    if (sum(grepl("/",getIndNames(gin)))==length(getIndNames(gin)))
    {
        longpops=getStrata(gin)
        print("rightslash found in Indnames")
    } else {
        longpops=sapply(strsplit(getIndNames(gin),"_"),function(x) x[1])
        print("rightslash absent in Indnames")
    }
#    print(longpops)
    names(longpops) = getIndNames(gin)
    orig.strata = getStrata(gin)
    setStrata(gin) = longpops
    popHet = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]
    names(popHet) = paste0(meta$pop,"Het")
    popPW=pairwiseTest(gin,nrep=0,stats="Fst_prime",quiet=T)
    popPair=sapply(popPW,function(x) x$result["Fst",1])
    names(popPair) = sapply(popPW,function(x) paste(names(x$strata.freq),collapse="_"))
    

    
    intro = meta[longpops,"intro"]
    names(intro) = getIndNames(gin)
    setStrata(gin) = intro
    intro_vs_nativeHet = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]
    names(intro_vs_nativeHet) = c("nativeHet","introducedHet")
    intro_vs_nativeOverall = overallTest(gin,nrep=0)$result[,1]
    names(intro_vs_nativeOverall) = paste0(names(intro_vs_nativeOverall),".intro")
    introPW=pairwiseTest(gin,nrep=0,stats="Fst_prime",quiet=T)
    intro_vs_nativePair=sapply(introPW,function(x) x$result["Fst",1])
    names(intro_vs_nativePair) = sapply(introPW,function(x) paste(names(x$strata.freq),collapse="_"))
    
    reg=meta[longpops,"source"]
    names(reg) = getIndNames(gin)
    setStrata(gin) = reg
    regHet = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]
    names(regHet) =  unique(reg)
    regOverall = overallTest(gin,nrep=0)$result[,1]
    names(regOverall) = paste0(names(regOverall),".reg")
    regPW=pairwiseTest(gin,nrep=0,stats="Fst_prime",quiet=T)
    regPair=sapply(regPW,function(x) x$result["Fst",1])
    names(regPair) = sapply(regPW,function(x) paste(names(x$strata.freq),collapse="_"))
    r=c(overallHet, ostruct,popHet,popPair,intro_vs_nativeHet, intro_vs_nativeOverall,
      intro_vs_nativePair,regHet,regOverall,regPair)
    r[order(names(r))]
}

###use the trad version for production
summary_stats_snp = summary_stats_trad


####
#### summary stats microsatellites
####


summary_stats.new = function(gin,meta)
{

    meta.orig=meta
    rownames(meta)=meta$longpop
    overallHet = mean(heterozygosity(gin)[,2])
    names(overallHet)="overallHet"

    if (sum(grepl("/",getIndNames(gin)))==length(getIndNames(gin)))
    {
        longpops=getStrata(gin)
    } else {
        longpops=sapply(strsplit(getIndNames(gin),"_"),function(x) x[3])
    }
    

names(longpops) = getIndNames(gin)
orig.strata = getStrata(gin)
setStrata(gin) = longpops
popHet = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]

names(popHet) = paste0(meta$pop,"Het")

popOverall = overallTest(gin,nrep=0)$result[,1]
names(popOverall) = paste0(names(popOverall),".pop")
#popPW=pairwiseTest(gin,nrep=0,stats="Fst_prime",quiet=T)
#popPair=sapply(popPW,function(x) x$result["Fst",1])
#names(popPair) = sapply(popPW,function(x) paste(names(x$strata.freq),collapse="_"))
    
intro = meta[longpops,"intro"]
names(intro) = getIndNames(gin)
setStrata(gin) = intro
intro_vs_nativeHet = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]
names(intro_vs_nativeHet) = c("nativeHet","introducedHet")
intro_vs_nativeOverall = overallTest(gin,nrep=0)$result[,1]
names(intro_vs_nativeOverall) = paste0(names(intro_vs_nativeOverall),".intro")
introPW=pairwiseTest(gin,nrep=0,stats="Fst_prime",quiet=T)
intro_vs_nativePair=sapply(introPW,function(x) x$result["Fst",1])
names(intro_vs_nativePair) = sapply(introPW,function(x) paste(names(x$strata.freq),collapse="_"))

reg=meta[longpops,"gigas_source"]
names(reg) = getIndNames(gin)
setStrata(gin) = reg
regHet = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]
names(regHet) =  unique(reg)
regM = with(mRatio(gin,by.strata=T),aggregate(cbind(m.ratio=m.ratio),by=list(stratum=stratum),mean,na.rm=T))[,2]
names(regM) =  unique(reg)
regOverall = overallTest(gin,nrep=0)$result[,1]
names(regOverall) = paste0(names(regOverall),".reg")
regPW=pairwiseTest(gin,nrep=0,stats="Fst_prime",quiet=T)
regPair=sapply(regPW,function(x) x$result["Fst",1])
names(regPair) = sapply(regPW,function(x) paste(names(x$strata.freq),collapse="_"))



reg=meta[longpops,"reg"]
names(reg) = getIndNames(gin)
setStrata(gin) = reg
reg2Het = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]
reg2M = with(mRatio(gin,by.strata=T),aggregate(cbind(m.ratio=m.ratio),by=list(stratum=stratum),mean,na.rm=T))[,2]
reg2Overall = overallTest(gin,nrep=0)$result[,1]
names(reg2Overall) = paste0(names(reg2Overall),".reg2")
reg2PW=pairwiseTest(gin,nrep=0,stats="Fst_prime",quiet=T)
reg2Pair=sapply(reg2PW,function(x) x$result["Fst",1])
names(reg2Pair) = sapply(reg2PW,function(x) paste(names(x$strata.freq),collapse="_"))


c(overallHet,popHet,intro_vs_nativeHet,regHet,reg2Het,
    popOverall,intro_vs_nativeOverall,regOverall,reg2Overall,
  intro_vs_nativePair,regPair,reg2Pair
  )
}


summary_stats.orig = function(gin,meta)
{

    meta.orig=meta
    rownames(meta)=meta$longpop
    overallHet = mean(heterozygosity(gin)[,2])
    names(overallHet)="overallHet"

    if (sum(grepl("/",getIndNames(gin)))==length(getIndNames(gin)))
    {
        longpops=getStrata(gin)
    } else {
        longpops=sapply(strsplit(getIndNames(gin),"_"),function(x) x[3])
    }
    

names(longpops) = getIndNames(gin)
orig.strata = getStrata(gin)
setStrata(gin) = longpops
popHet = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]
popM = with(mRatio(gin,by.strata=T),aggregate(cbind(m.ratio=m.ratio),by=list(stratum=stratum),mean,na.rm=T))[,2]
names(popHet) = paste0(meta$pop,"Het")
names(popM) = paste0(meta$pop,"M")
popOverall = overallTest(gin,nrep=0)$result[,1]
names(popOverall) = paste0(names(popOverall),".pop")
#popPW=pairwiseTest(gin,nrep=0,stats="Fst_prime",quiet=T)
#popPair=sapply(popPW,function(x) x$result["Fst",1])
#names(popPair) = sapply(popPW,function(x) paste(names(x$strata.freq),collapse="_"))
    
intro = meta[longpops,"intro"]
names(intro) = getIndNames(gin)
setStrata(gin) = intro
intro_vs_nativeHet = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]
names(intro_vs_nativeHet) = c("nativeHet","introducedHet")
intro_vs_nativeM = with(mRatio(gin,by.strata=T),aggregate(cbind(m.ratio=m.ratio),by=list(stratum=stratum),mean,na.rm=T))[,2]
names(intro_vs_nativeHet) = c("nativeM","introducedM")
intro_vs_nativeOverall = overallTest(gin,nrep=0)$result[,1]
names(intro_vs_nativeOverall) = paste0(names(intro_vs_nativeOverall),".intro")
introPW=pairwiseTest(gin,nrep=0,stats="Fst_prime",quiet=T)
intro_vs_nativePair=sapply(introPW,function(x) x$result["Fst",1])
names(intro_vs_nativePair) = sapply(introPW,function(x) paste(names(x$strata.freq),collapse="_"))

reg=meta[longpops,"gigas_source"]
names(reg) = getIndNames(gin)
setStrata(gin) = reg
regHet = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]
names(regHet) =  unique(reg)
regM = with(mRatio(gin,by.strata=T),aggregate(cbind(m.ratio=m.ratio),by=list(stratum=stratum),mean,na.rm=T))[,2]
names(regM) =  unique(reg)
regOverall = overallTest(gin,nrep=0)$result[,1]
names(regOverall) = paste0(names(regOverall),".reg")
regPW=pairwiseTest(gin,nrep=0,stats="Fst_prime",quiet=T)
regPair=sapply(regPW,function(x) x$result["Fst",1])
names(regPair) = sapply(regPW,function(x) paste(names(x$strata.freq),collapse="_"))



reg=meta[longpops,"reg"]
names(reg) = getIndNames(gin)
setStrata(gin) = reg
reg2Het = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]
reg2M = with(mRatio(gin,by.strata=T),aggregate(cbind(m.ratio=m.ratio),by=list(stratum=stratum),mean,na.rm=T))[,2]
reg2Overall = overallTest(gin,nrep=0)$result[,1]
names(reg2Overall) = paste0(names(reg2Overall),".reg2")
reg2PW=pairwiseTest(gin,nrep=0,stats="Fst_prime",quiet=T)
reg2Pair=sapply(reg2PW,function(x) x$result["Fst",1])
names(reg2Pair) = sapply(reg2PW,function(x) paste(names(x$strata.freq),collapse="_"))


c(overallHet,popHet,intro_vs_nativeHet,regHet,reg2Het,
  popM,intro_vs_nativeM,regM,reg2M,
  popOverall,intro_vs_nativeOverall,regOverall,reg2Overall,
  intro_vs_nativePair,regPair,reg2Pair
  )
}


#' Summary Statistics for Microsatellites
#'
#' This function computes various summary statistics related to microsatellite data.
#'
#' @param gin Genetic input data in strataG gtypes format
#' @param meta Metadata related to the `gin` input.
#'
#' @details 
#' The function computes multiple summary statistics such as overall heterozygosity, 
#' population-specific heterozygosity, pairwise population statistics, and so on.
#' The method to extract population information from `gin` is also determined within the function based on 
#' the presence of a specific pattern in the individual names.
#'
#' @return A vector containing various summary statistics. The output is ordered to ensure comparability 
#' between original and simulated statistics.
#'
#' @seealso 
#' \code{\link{heterozygosity}}, \code{\link{overallTest}}, \code{\link{pairwiseTest}}, \code{\link{mRatio}}
#'

summary_stats_microsatellite = function(gin,meta)
{

    meta.orig=meta
    meta$longpop=meta$pop
    rownames(meta)=meta$longpop
    meta$intro = ifelse(meta$source %in% c("PNW","Cali","EU"),T,F)
    overallHet = mean(heterozygosity(gin)[,2])
    names(overallHet)="overallHet"

    if (sum(grepl("/",getIndNames(gin)))==length(getIndNames(gin)))
    {
        longpops=getStrata(gin)
    } else {
        longpops=gsub("[0-9]+","",getIndNames(gin))#sapply(strsplit(getIndNames(gin),"_"),function(x) x[3])
    }
    

names(longpops) = getIndNames(gin)
orig.strata = getStrata(gin)
setStrata(gin) = longpops
popHet = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]

names(popHet) = paste0(meta$pop,"Het")

popOverall = overallTest(gin,nrep=0)$result[,1]
names(popOverall) = paste0(names(popOverall),".pop")

    d=as.matrix(adegenet::dist.genpop(adegenet::genind2genpop(strataG::gtypes2genind(gin),quiet=T)))
    pwdf = unique(t(apply(expand.grid(col=colnames(d),row=rownames(d)),1,sort)))
    pwdf = pwdf[pwdf[,1]!=pwdf[,2],]
    colnames(pwdf)=c("col","row")
    popPair=sapply(1:nrow(pwdf),function(x){r=d[pwdf[x,1],pwdf[x,2]];names(r)=paste0(pwdf[x,1],'_',pwdf[x,2]);r})
    
intro = meta[longpops,"intro"]
names(intro) = getIndNames(gin)
setStrata(gin) = intro
intro_vs_nativeHet = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]
names(intro_vs_nativeHet) = c("nativeHet","introducedHet")
intro_vs_nativeOverall = overallTest(gin,nrep=0)$result[,1]
names(intro_vs_nativeOverall) = paste0(names(intro_vs_nativeOverall),".intro")
introPW=pairwiseTest(gin,nrep=0,stats="Fst_prime",quiet=T)
intro_vs_nativePair=sapply(introPW,function(x) x$result["Fst",1])
names(intro_vs_nativePair) = sapply(introPW,function(x) paste(names(x$strata.freq),collapse="_"))

reg=meta[longpops,"source"]
names(reg) = getIndNames(gin)
setStrata(gin) = reg
regHet = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]
names(regHet) =  unique(reg)
regM = with(mRatio(gin,by.strata=T),aggregate(cbind(m.ratio=m.ratio),by=list(stratum=stratum),mean,na.rm=T))[,2]
names(regM) =  unique(reg)
regOverall = overallTest(gin,nrep=0)$result[,1]
names(regOverall) = paste0(names(regOverall),".reg")
regPW=pairwiseTest(gin,nrep=0,stats="Fst_prime",quiet=T)
regPair=sapply(regPW,function(x) x$result["Fst",1])
names(regPair) = sapply(regPW,function(x) paste(names(x$strata.freq),collapse="_"))

    rv=c(overallHet,popHet,intro_vs_nativeHet,intro_vs_nativeOverall,
         intro_vs_nativePair,regHet,regM,regOverall,regPair,popPair)
    rv=rv[order(names(rv))]  #make sure that the original and simulated summary stats are directly comparable
    rv
}
