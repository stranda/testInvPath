library(testInvPath)
library(abc)
library(abcrf)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)

### Check if an argument is provided, and if it's numeric
if (length(args) > 0 && !is.na(as.numeric(args[1]))) {
    cores <- as.numeric(args[1])
} else {
    cores <- 1
}


diagnosticPlot = FALSE
crossVal = FALSE  #set to true to run (really slow) cross-validations
source("datafiles.R")



##read in the raw summary stats, process and make a reference table

genofile=paste(abspath,genofile,sep='/')
mname=paste(abspath,mname,sep='/')
utfn=paste0(abspath,'/data/untransRef.RDS')

if (!file.exists(utfn))
{
    untrans.ref=make_ref_table(genofile=genofile,
                               mname=mname,
                               intros=intros,sources=sources,
                               csvpath=paste0(abspath,'/gathered'),
                               dataType=dataType,
                               popPairwise=popPairwise)
    saveRDS(file=utfn,untrans.ref)
} else untrans.ref = readRDS(file=utfn)

##########################
##
## pca transform the reference table and observations
##
pcafn=paste0(abspath,'/data/pcaRef.RDS')
if (!file.exists(pcafn))
    {
        pca.ref = make_pca_ref_table(untrans.ref,prop.variation=0.99)
        saveRDS(file=pcafn,pca.ref)
    } else pca.ref=readRDS(file=pcafn)


###measures of model quality
###first plot the observed in the rest of PCA space
pdf(file=paste0(abspath,"/figs/pca.pdf"))
plotPCA(pca.ref,prop_points=0.2)
dev.off()


###examine the differences between summary stat means and the observed data
pdf(file=paste0(abspath,"/figs/obsVsSimsSumStats.pdf"))
obsModelDist(untrans.ref,plot=T) ->tmp
dev.off()

################# estimate posterior probs for models
##### first the scaled summary stats
scssfn=paste0(abspath,'/data/modsel_untrans.RDS')
if (!file.exists(scssfn))
    {
        modsel_untr = ModelPostProbs(untrans.ref,method=c("mnlogistic"),
                                     cores=cores,modelColumn="compositeModel")
        saveRDS(file=scssfn,modsel_untr)
    } else modsel_untr=readRDS(file=scssfn)

###cross validation
cvfn=paste0(abspath,'/data/cv4postpr_untrans.RDS')
if (crossVal)
    if (!file.exists(cvfn))
    {
        cv_untr = CVModelPostProbs(untrans.ref,method=c("mnlogistic"),
                                     cores=cores,modelColumn="compositeModel")
        saveRDS(file=cvfn,cv_untr)
    } else cv_untr=readRDS(file=cvfn)


##
## model selection using the pca transformed reference table
##

pcassfn=paste0(abspath,'/data/modsel_pca.RDS')
if (!file.exists(pcassfn))
    {
        modsel_pca = ModelPostProbs(pca.ref,method=c("mnlogistic","neuralnet"),
                                     cores=cores,modelColumn="compositeModel")
        saveRDS(file=pcassfn,modsel_pca)
    } else modsel_pca=readRDS(file=pcassfn)

###cross validation
cvpfn=paste0(abspath,'/data/cv4postpr_pca.RDS')
if (crossVal)
    if (!file.exists(cvpfn))
    {
        cv_pca = CVModelPostProbs(pca.ref,method=c("mnlogistic"),
                                     cores=cores,modelColumn="compositeModel")
        saveRDS(file=cvpfn,cv_untr)
    } else cv_pca=readRDS(file=cvpfn)

###
probdf = data.frame(do.call(rbind,lapply(modsel_pca, function(x)
{
    y=summary(x,print=F,rejection=F)$Prob
})))

probdf$method=sapply(strsplit(rownames(probdf),"_"),function(x)x[1])
probdf$tol=as.numeric(sapply(strsplit(rownames(probdf),"_"),function(x)x[2]))

pdf(file=paste0(abspath,"/figs/posterior_pca.pdf"))
barplot(as.matrix(probdf[probdf$method=="neuralnet",1:8]),beside=T,names=gsub("X","",names(probdf)[1:8]),ylab="posterior probability of model",xlab="model (1-2==shipping, 3-4==gigas)",main=species)
legend(x=1,y=0.9*max(c(unlist(probdf[probdf$method=="neuralnet",1:8]))),legend=probdf[probdf$method=="neuralnet","tol"],title="Tol",fill=grey.colors(6)[1:6])
dev.off()

ship_v_gigas = cbind(ship_pr=rowSums(probdf[probdf$method=="neuralnet",1:4]),gigas_pr=rowSums(probdf[probdf$method=="neuralnet",5:8]))
print(ship_v_gigas)
pdf(file=paste0(abspath,"/figs/bayes_factor_pca.pdf"))
barplot(as.matrix(ship_v_gigas),beside=T,main=species)
dev.off()

probdf$species=species
write.table(file=paste0(abspath,"/data/probdf.csv"),sep=",",row.names=F,probdf)

##
## build a random forest of the summary stats and train to predict the model indices
## run it on PCA transformed vars
##
rffn=paste0(abspath,'/data/modelRF_pca.RDS')
if (!file.exists(rffn))
    {
        modelRF_pca = modelRF(pca.ref,cores=cores,modelColumn="compositeModel")
        saveRDS(file=rffn,modelRF_pca)
    } else modelRF_pca=readRDS(file=rffn)


rfpfn=paste0(abspath,'/data/modelRF_pca_pred.RDS')
if (!file.exists(rfpfn))
    {
        modelRF_pca_pred = predRFmodel(modelRF_pca,pca.ref,cores=cores,modelColumn="compositeModel")
        saveRDS(file=rfpfn,modelRF_pca_pred)
    } else modelRF_pca_pred=readRDS(file=rfpfn)


if (diagnosticPlot)
{
pdf(file=paste0(abspath,"/figs/modelRF_pca.pdf"))
plotRFmodel(modelRF_pca,pca.ref)
dev.off()

pdf(file=paste0(abspath,"/figs/modelRF_pca_err.pdf"))
modelRF_pca_err = errRFmodel(modelRF_pca,pca.ref,cores=cores)
dev.off()
}

##
## build a random forest of the summary stats and train to predict the model indices
## run it on scaled but untransformed vars
##
rffn=paste0(abspath,'/data/modelRF.RDS')
if (!file.exists(rffn))
    {
        modelRF = modelRF(untrans.ref,cores=cores,modelColumn="compositeModel")
        saveRDS(file=rffn,modelRF)
    } else modelRF=readRDS(file=rffn)


rfpfn=paste0(abspath,'/data/modelRF_pred.RDS')
if (!file.exists(rfpfn))
    {
        modelRF_pred = predRFmodel(modelRF,untrans.ref,cores=cores,modelColumn="compositeModel")
        saveRDS(file=rfpfn,modelRF_pred)
    } else modelRF_pred=readRDS(file=rfpfn)


if (diagnosticPlot)
    {
###RF model diagnostics

        pdf(file=paste0(abspath,"/figs/modelRF.pdf"))
        par(ask=FALSE)
        plotRFmodel(modelRF,untrans.ref)
        dev.off()
        
        pdf(file=paste0(abspath,"/figs/modelRF_err.pdf"))
        par(ask=FALSE)
        modelRF_err = errRFmodel(modelRF,untrans.ref,cores=cores)
        dev.off()
    }



################ parameter estimates

paramfn=paste0(abspath,'/data/params_untrans.RDS')
if (!file.exists(paramfn))
{
    post_params_untrans = paramPosteriors(untrans.ref,method=c("loclinear","ridge","neuralnet"),cores=cores,tol=c(0.1,0.01))
    saveRDS(file=paramfn,post_params_untrans)
} else post_params_untrans=readRDS(file=paramfn)




paramfn=paste0(abspath,'/data/params_pca.RDS')
if (!file.exists(paramfn))
{
    post_params_pca = paramPosteriors(pca.ref,method=c("loclinear","ridge","neuralnet"),cores=cores,tol=c(0.1,0.05,0.01))
    saveRDS(file=paramfn,post_params_pca)
} else post_params_pca=readRDS(file=paramfn)

##Plot PCA-based parameter estimates

for (n in names(post_params_pca))
     {
         pdf(file=paste0(abspath,"/figs/",n,"_param_est_pca.pdf"))
         plot(post_params_pca[[n]],pca.ref$params[,1:14],ask=F)
         dev.off()
     }


##Plot untransformed summary stats-based parameter estimates

for (n in names(post_params_pca))
     {
         pdf(file=paste0(abspath,"/figs/",n,"_param_est.pdf"))
         plot(post_params_untrans[[n]],untrans.ref$params[,1:14],ask=F)
         dev.off()
     }


