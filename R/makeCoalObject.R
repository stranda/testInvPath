#' Create a single object ready for coalescent simulation based on pulls from prior distributions
#'
#' @param params these are the parameters produced by setupReps()
#' @param priors pulls from prior distribution to parameterize this
#'     rep (from getPriors())
#' @param demeN the total number of demes to be simulated in the
#'     native (source range)
#' @param demeI the total number of demes to be simulated in the
#'     introduced (source range)
#' @description this function currently uses strataG to parameterize a
#'     fastsimcoal model that depends on the structure determined in
#'     params (primarily source population, native topology,
#'     introduced populations).  the priors determine choices of
#'     timing, actual introduction model and rates of mutation and
#'     migration.  demeN and demeI try to simulate ghost populations
#'     that are ultimately not sampled.  If the values for these
#'     variables are less that the number of sampled demes in each
#'     region then just the sampled populations are simulated.
#'     Otherwise additional unsampled populations are simulated up to
#'     the number of demeN or demeI for each of the native or
#'     introduced regions, respectively
#' @return a fscWrite() object that can be used in strataG's fscRun() function
#' @export

createCoalObject <- function(params,priors,demeN=0,demeI=0)
{
    meta=params$meta
    meta=meta[order(meta$longpop),]
    meta$idnum=(1:nrow(meta))-1

### inset ghost populations into the meta object. If demeI or demeN are zero then no ghosts created in intro or native, respectively
    if (max(c(demeN,demeI))>0)
    {
        regtbl = meta%>%group_by(intro,source)%>%summarize(n=n()) #this dataframe contains the number of sampled pops in each region in the metafile
        madd=do.call(rbind,lapply(unique(meta$source),function(s)  #this assumes no population is both source and introduced (not modeling hopping)
        {
            print(s)
            if (regtbl$intro[regtbl$source==s]) #these are introduced pops
            {
                iNeeded = demeI-regtbl$n[regtbl$source==s]
                if (iNeeded>0)
                {
                    m=meta[rep(which(meta$source==s)[1],iNeeded),]
                }
            } else { #these would be the native pops
                nNeeded = demeN-regtbl$n[regtbl$source==s]
                if (nNeeded>0)
                {
                    m=meta[rep(which(meta$source==s)[1],nNeeded),]
                } 
            }
            if (exists("m")) #this can happen if there are more sampled pops in a region than demes we ask for
                {
                    m$longpop=rep(paste0(s,"_ghost"),nrow(m))
                    m$idnum=rep(NA,nrow(m))
                    m
                } else {
                    NULL
                }
        }))
        madd$longpop=paste0(madd$longpop,1:nrow(madd))
        madd$pop=sprintf("%03i",1:nrow(madd))
        madd$idnum=max(meta$idnum)+(1:nrow(madd))
        meta=rbind(meta,madd)
    }
    
    itbl = params$introTbls[[priors$introModel]]
#    print(dim(itbl))
    itbl = apply(itbl,2,function(x){x/sum(x)})
#    print(dim(itbl))
    itbl[is.na(itbl)]=0
    plmult = ifelse((params$dataType=="sequence")&(params$ploidy==2),2,1)
    poptbl = as.table(params$samples[,2]*plmult)
    names(poptbl)=params$samples[,1]
    if (length(poptbl)!=sum(!grepl("ghost",meta$longpop)))
        print("there are pops in meta that are not ghosts and also have no sample")

    ##alter poptbl to reflect the ghosts created above
     ghosts=rep(0,sum(grepl("ghost",meta$longpop)))
    names(ghosts)=meta$longpop[grepl("ghost",meta$longpop)]
    ghosts=as.table(ghosts)
    poptbl=c(poptbl,ghosts)

    demes = strataG::fscSettingsDemes(do.call(rbind,lapply(names(poptbl),function(p)
    {
        strataG::fscDeme(deme.size=ifelse(meta$intro[meta$longpop==p],priors$introNe,priors$nativeNe),
                sample.size=poptbl[p],
                growth=ifelse(meta$intro[meta$longpop==p],priors$IntroGrow,0))
    }
    )),ploidy=ifelse(params$dataType=="sequence",1,2))
    demes$deme.name = names(poptbl)

    
#### starting to sort out events
### make two lists, one has the sampled deme ids of the intropops and the other, native pops

#    print(itbl)
#    print(colnames(itbl))
    ipops = lapply(colnames(itbl),function(r)
    {
        meta$idnum[meta$source==r]
    }
    )
    names(ipops)=colnames(itbl)
    npops = lapply(rownames(itbl),function(r)
    {
        meta$idnum[meta$source==r]
    })
    names(npops)=rownames(itbl)
#print("printing ipops")
#    print(ipops)
#print("printing npops")
#    print(npops)
    
###events in the native range within regions
    nEventsWInReg = data.frame(do.call(rbind,lapply(names(npops),function(r)
    {
        dn = npops[[r]]
        if (length(dn)>1)
        {
            m=matrix(c(priors$t0,NA,dn[1],1,1,0,1),
                     byrow=T,ncol=7,nrow=length(dn)-1)
            m[,2]=dn[-1]
            m2=m #try to set deme sizes to zero after coalescence
            m2[,3]=m2[,2]
            m2[,4:6]=0
            m2[,7]=2 #mig mat 1 is all zeros
            
            data.frame(rbind(m,m2))
        }
    })))
###events in the native range among regions
#    print("this far 1")
#    print("about to print npops")
#    print(npops)
    nEventsAmongReg=data.frame(do.call(rbind,lapply(names(npops),function(r)
    {
#        print(r)
        src = npops[[r]][1] #each region has already coalesced into the first deme, now these coalesce
#        print(npops)
#        print(src)
        nrt=params$native_range_topology
#        print(nrt)
#        print(class(nrt))
        if (max(nrt[,r])>0)
        {
            sink=npops[[names(which(nrt[,r]>0))]][1]
#            print(paste(src,sink))
            timeflag=nrt[names(which(nrt[,r]>0)),r]
            ret1=c(c(priors$t1,priors$t2,priors$t3)[timeflag],  ##this chooses the correct time for the particular coalescence...
                   src, sink, 1,
                   ifelse(timeflag!=3,1,priors$ancestralNe), #new size (in past) 1 unless final coalescence
                   0,
                   2)
            ret2=ret1   ##ret2 is supposed to set abandonded demes to size 0
 #           ret2[1]=ret1[1]+1 #add a gen before extinction
            ret2[4:6]=0
            ret2[7]=2
            ret2[3]=ret2[2] 
            rbind(ret1,ret2)
        }
    })))
    nEventsAmongReg = nEventsAmongReg[order(nEventsAmongReg[,1]),]
#    print("this far")
#    print(nEventsAmongReg)
    
    iEventsIntro = do.call(rbind,lapply(names(ipops),function(r)
    {
        sources=ipops[[r]]
        ivec=itbl[,r]
        ivec=ivec[ivec>0]

        if (length(ivec)>1)
            if (!priors$SimulIntro)
            {
                ivec=sample(ivec,1,prob=ivec)
                ivec[1]=1
            }
        
        rsrc=sapply(names(npops),function(n){if (length(npops[[n]])>1) sample(npops[[n]],1) else npops[[n]]})#which pop in native range to coalesce into
        
#        if (meta$source[meta$idnum%in%rsrc]%in%names(npops)) print ("introduced coalescing into introduced")
        names(rsrc)=names(npops)
#print(r)
        ret1=do.call(rbind,lapply(names(ivec),function(sr)
        {
            m=matrix(c(NA),ncol=7,nrow=length(sources))
            m[,1]=priors$tIntro-1 #assign times
            m[,2]=sources
            m[,3]=rsrc[sr]
            ##the next line assigns the proportion of lineages to ivec.  If the last event, all lineages coalesce
            if (which(names(ivec)==sr)!=length(ivec)) m[,4]=ivec[sr] else m[,4]=1
            m[,5]=1 #new size
            m[,6]=0
            m[,7]=0 #mig matrix
            m=data.frame(m)
            m[m[,2]!=m[1,2],3]=sources[1] #locally, everybody coalesces into sources[1]
            m[m[,2]==m[1,2],1]=priors$tIntro #this is when sources[1] coalesces into native
            m
        }))

        ret2=ret1
        ret2[,4:6]=0
        ret2[,7]=1
        ret2[,3]=ret2[,2]
        rbind(ret1,unique(ret2))
    
    }))

    edf = data.frame(rbind(iEventsIntro, nEventsWInReg,  nEventsAmongReg))

    edf = edf[order(edf[,1],edf[,2],edf[,2]==edf[,3],edf[,3]),]

    events = strataG::fscSettingsEvents(do.call(rbind,lapply(1:nrow(edf),function(r)
    {
        strataG::fscEvent(round(edf[r,1],2),edf[r,2],edf[r,3],edf[r,4],edf[r,5],edf[r,6],edf[r,7])
    }
    )))
#    events$event.time = ifelse(events$source==events$sink,events$event.time+0.01,events$event.time)
#print(events)    
    basemat = matrix(0,nrow(meta),nrow(meta))
    mig.mat = basemat
    for (i in 1:nrow(mig.mat))
        for (j in 1:nrow(mig.mat))
            if (i!=j)
            {
                src.intro=meta$intro[i]
                src.reg=meta$source[i]
                snk.intro=meta$intro[j]
                snk.reg=meta$source[j]
                if (src.reg==snk.reg)
                    if (snk.intro & src.intro)
                        mig.mat[i,j]=priors$introM else mig.mat[i,j]=priors$nativeM
            }
    mig.mat.native = basemat
    for (i in 1:nrow(mig.mat.native))
        for (j in 1:nrow(mig.mat.native))
            if (i!=j)
            {
                src.intro=meta$intro[i]
                src.reg=meta$source[i]
                snk.intro=meta$intro[j]
                snk.reg=meta$source[j]
                if (src.reg==snk.reg)
                    if (snk.intro & src.intro)
                        mig.mat.native[i,j]=0 else mig.mat.native[i,j]=priors$nativeM
            }
    
    migration = strataG::fscSettingsMigration(mig.mat,mig.mat.native,basemat)

    genetics <- switch(params$dataType,
                 "sequence" = {
                     ## simulate DNA sequence
                     ## set the length of the sequence to use
                     seqlen=ncol(params$gin@sequences@dna$gene.1)
                     strataG::fscSettingsGenetics(
                                  strataG::fscBlock_dna(sequence.length=seqlen,mut.rate=priors$mut_rate,transition.rate=0.33),
                                  strataG::fscBlock_dna(sequence.length=seqlen,mut.rate=priors$mut_rate,transition.rate=0.33),
                                  num.chrom=NULL)
                 },
                 "snp" = {
                   # simulate snps
                     strataG::fscSettingsGenetics(
                                  strataG::fscBlock_snp(sequence.length=1,
                                                        mut.rate=priors$mut_rate),
                                  num.chrom=1000)
                 },
                 "microsatellite" = {
                     strataG::fscSettingsGenetics(
                         strataG::fscBlock_microsat(1,priors$mut_rate,0,chromosome=1)
                     )
                 },
                 {
                   # Default action in case no match is found
                   stop("params$dataType must be 'sequence', 'snp', or 'microsatellite'")
                 }
                 )
    if (params$dataType=="microsatellite") #we are going to set the correct number of genetic loci now
    {
        nl=getNumLoci(params$gin)
        genetics=genetics[rep(1,nl),]
        genetics$chromosome=c(1:nl)
    }
    

##now assemble all of the components into a single set of parameters and return

    strataG::fscWrite(demes,genetics,migration,events,use.wd=F)
}

