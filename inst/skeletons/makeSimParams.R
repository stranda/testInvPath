makeSimParams <- function()
{

equalSourceTbl = gigasTblML
equalSourceTbl[!is.na(equalSourceTbl)]=1

introTbls = list(shipTbl30=shipTbl30,shipTbl60=shipTbl60,gigasTblAdmix=gigasTblAdmix,gigasTblML=gigasTblML,equalSourceTbl=equalSourceTbl)



##############################################################
meta1 = read.csv(metaname1)
meta1$longpop=meta1$pop
gin=readRDS(genofile)
setOther(gin,"popnames")=getStrata(gin)
gin.orig=gin

###undaria specific
meta1$source[meta1$pop %in% c("Pop38","Pop39","Pop40","Pop41")]="Cali"
meta1$source[meta1$pop %in% c("Pop70","Pop71")]="EU"
###

rmpops="empty"
##added for undaria, might not need
rmpops = meta1$pop[grep("China|nonSource",meta1$source)]
###
keepvec = !(getStrata(gin)%in% rmpops)
gin=gin[,,keepvec]#remove non-source
meta1=meta1[!(meta1$pop %in% rmpops),]

###keep a lot of different possible strata with the gin object
###indmeta is just the meta1 df mapped onto the gridID currently used as strata
indStrat=data.frame(ind=getIndNames(gin),pop=getStrata(gin))
indmeta=merge(indStrat , meta1 , all.x=T)
rownames(indmeta) = getIndNames(gin)
setOther(gin,"indmeta")=indmeta
    

meta = meta1
meta$idnum = 0:(nrow(meta)-1)


pop=getOther(gin,"popnames")
pop=pop[pop %in% meta$pop]
demedf = data.frame(longpop=unique(pop))

samples = as.data.frame(table(pop))
names(samples) = c("longpop","sample.size")
meta$intro = meta$reg != "1_Asia"
meta = meta[,c("longpop","intro","source","idnum","pop")]

#meta = meta[gsub(" ","",meta$pop) %in% gsub(" ","",unique(getStrata(gin))),]


##the native range topology maps the times t1-t3 to the 4 regions in japan
##will use this in the generate coalescence model functions to determine within
##native range coalescence.  This works by going down columns of the matrix.  Hokkaido
##coalesces into honshu at time 1 and honshu into tokyo at time 3 (at least for the
##original gvermSNP file.  if this is downstream, regions likely differ

                                #k,h,s,t,hk    
native_range_topology = matrix(c(0,0,0,0,0,  #kag
                                 0,0,0,0,1,  #hon
                                 1,0,0,0,0,  #sea
                                 0,3,2,0,0,  #tok
                                 0,0,0,0,0), #hok
                                 byrow=T,nrow=5)
                               
colnames(native_range_topology)=c("kag","hon","sea","tok","hok")
rownames(native_range_topology)=colnames(native_range_topology)

list(gin=gin,species=species,meta=meta,introTbls=introTbls,native_range_topology=native_range_topology,demeN=demeN,demeI=demeI,samples=samples)
}



