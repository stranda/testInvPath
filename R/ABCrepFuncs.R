#' Run ABC Prior Replicate Analysis with Diagnostic and Model Selection
#'
#' This function performs Approximate Bayesian Computation (ABC) using prior probabilities for different scenarios,
#' transforms data using PCA, runs model selection, performs cross-validation if enabled, and builds Random Forest models.
#' It also generates various diagnostic plots to evaluate the performance of the models.
#'
#' @param probShip Numeric, the prior probability of shipping (default is 0.5).
#' @param priortag String, a tag for naming results from this prior analysis.
#' @param priorrep Integer, the replicate number.
#' @param datafn String, the filename of the R script that loads necessary data.
#' @param diagnosticPlot Logical, flag to generate diagnostic plots (default is FALSE).
#' @param diagnosticPlotRF Logical, flag to generate Random Forest diagnostic plots (default is FALSE).
#' @param crossVal Logical, whether to perform cross-validation (default is FALSE; note: this can be very slow).
#' @param rfsample Integer, the number of samples for the Random Forest (default is 150000).
#' @param cores Integer, the number of cores to use for parallel processing (default is 1).
#' @param verb Logical, whether to print verbose messages (default is TRUE).
#' @return A list of paths to files generated during the process and model results.
#' @examples
#' runABCpriorRep(probShip = 0.7, priortag = "testPrior", priorrep = 1, datafn = "data_loader.R",
#'                diagnosticPlot = TRUE, diagnosticPlotRF = TRUE, crossVal = FALSE)
#' @export
runABCpriorRep = function(probShip=0.5,     #prior prob of shipping 
                          priortag="prior", #tag for naming results from this prior
                          priorrep=1,
                          datafn="datafiles.R",
                          diagnosticPlot = FALSE,
                          diagnosticPlotRF = FALSE,
                          crossVal = FALSE, #set to true to run (really slow) cross-validations
                          rfsample=150000,
                          cores=1,
                          verb=TRUE)

{
### initialize
    source(datafn)

    ##setup file names
    prioradd = paste0(priortag,"_",priorrep)
    
    utfn=paste0(abspath,'/data/untransRef.RDS') #untrans ref file
    pcafn=paste0(abspath,'/data/pcaRef.RDS') #pca ref file
    pcassfn=paste0(abspath,'/data/modsel_pca_',prioradd,'.RDS') #results of pca model selection
    cvpfn=paste0(abspath,'/data/cv4postpr_pca_',prioradd,'.RDS') #pcr cv for modsel
    rffn=paste0(abspath,'/data/modelRF_pca_',prioradd,'.RDS') #random forest on pca transformed data
    rfpfn=paste0(abspath,'/data/modelRF_pca_pred_',prioradd,'.RDS') ##pca model preds

    rffnut=paste0(abspath,'/data/modelRF_',prioradd,'.RDS')
    rfpfnut=paste0(abspath,'/data/modelRF_pred_',prioradd,'.RDS')

    probdffn=paste0(abspath,'/data/probdf_',prioradd,'.csv') #write the probdf data to the disk (this is for standard abc)
    rfsumfn=paste0(abspath,'/data/rfsum_',prioradd,'.csv') #write out posterir stuff for RF
    rfsinkfn=paste0(abspath,'/data/RFmodelOutput_',prioradd,'.txt') ##sink() outout for rf
    pcaplotfn=paste0(abspath,"/figs/pca.pdf") ##make the plots for the pca versus observed (might not change with priors?)
    post_pca_fn=paste0(abspath,'/figs/posterior_pca_',prioradd,'.pdf') #plots of the posteriors from trad ABC
    bayesPlot_nn=paste0(abspath,'/figs/bayes_factor_pca_nn_',prioradd,'.pdf') ##bayes factors neural net
    bayesPlot_mnlog=paste0(abspath,'/figs/bayes_factor_pca_mnlog_',prioradd,'.pdf') ##bayes factors
    
    ##read in the raw summary stats, process and make a reference table
    genofile=paste(abspath,genofile,sep='/')
    mname=paste(abspath,mname,sep='/')


    ### start calcing!!

print(paste("untransfile name",utfn))
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
    print (paste(date(),"setup untrans ref"))

##########################
##
## pca transform the reference table and observations
    ##
print(paste("pca ref 
file name",pcafn))
    
if (!file.exists(pcafn))
    {
        pca.ref = make_pca_ref_table(untrans.ref,prop.variation=0.99)
        saveRDS(file=pcafn,pca.ref)
    } else pca.ref=readRDS(file=pcafn)
print (paste(date(),"setup pca ref"))

###measures of model quality
###first plot the observed in the rest of PCA space
pdf(file=pcaplotfn)
plotPCA(pca.ref,prop_points=0.2)
dev.off()

###
### make sure that the models adhere to the priors in the priorfile
###
### subset the models based on priors
###

print(table(untrans.ref$params$compositeModel))
print(table(untrans.ref$params$introModel))

rows2keep = c(sample(which(untrans.ref$params$introModel%in%c(1,2)),floor(probShip*rfsample)),
              sample(which(untrans.ref$params$introModel%in%c(3,4)),floor((1-probShip)*rfsample)))

print(length(rows2keep))
untrans.ref$params = untrans.ref$params[rows2keep,]
untrans.ref$ref = untrans.ref$ref[rows2keep,]

pca.ref$params = pca.ref$params[rows2keep,]
pca.ref$ref = pca.ref$ref[rows2keep,]

print(table(untrans.ref$params$compositeModel))
print(table(untrans.ref$params$introModel))

print(table(pca.ref$params$compositeModel))
print(table(pca.ref$params$introModel))

################# estimate posterior probs for models
##
## model selection using the pca transformed reference table
##
if (!file.exists(pcassfn))
    {
        modsel_pca = ModelPostProbs(pca.ref,method=c("mnlogistic","neuralnet"),
                                     cores=cores,modelColumn="compositeModel")
        saveRDS(file=pcassfn,modsel_pca)
    } else modsel_pca=readRDS(file=pcassfn)
print (paste(date(),"ran pca postprods"))

###cross validation

if (crossVal)
{
    if (!file.exists(cvpfn))
    {
        cv_pca = CVModelPostProbs(pca.ref,method=c("mnlogistic","neuralnet"),
                                  tol=c(0.05,0.01,0.005,0.001),
                                  cores=cores,modelColumn="compositeModel")
        saveRDS(file=cvpfn,cv_untr)
    } else cv_pca=readRDS(file=cvpfn)
print (paste(date(),"ran pca CV postprods"))
}
###

print(length(modsel_pca))
probdf = data.frame(do.call(rbind,lapply(modsel_pca, function(x)
{
    y=summary(x,print=F,rejection=F)$Prob
})))

print("made it 1 in abcprior function")
probdf$method=sapply(strsplit(rownames(probdf),"_"),function(x)x[1])
probdf$tol=as.numeric(sapply(strsplit(rownames(probdf),"_"),function(x)x[2]))

pdf(file=post_pca_fn)
barplot(as.matrix(probdf[probdf$method=="neuralnet",1:8]),beside=T,names=gsub("X","",names(probdf)[1:8]),ylab="posterior probability of model",xlab="model (1-2==shipping, 3-4==gigas)",main=paste(species,"NN"))
legend(x=1,y=0.9*max(c(unlist(probdf[probdf$method=="neuralnet",1:8]))),legend=probdf[probdf$method=="neuralnet","tol"],title="Tol",fill=grey.colors(6)[1:6])
barplot(as.matrix(probdf[probdf$method=="mnlogistic",1:8]),beside=T,names=gsub("X","",names(probdf)[1:8]),ylab="posterior probability of model",xlab="model (1-2==shipping, 3-4==gigas)",main=paste(species,"MnLog"))
legend(x=1,y=0.9*max(c(unlist(probdf[probdf$method=="mnlogistic",1:8]))),legend=probdf[probdf$method=="mnlogistic","tol"],title="Tol",fill=grey.colors(6)[1:6])
dev.off()

ship_v_gigas = cbind(ship_pr=rowSums(probdf[probdf$method=="neuralnet",1:4]),gigas_pr=rowSums(probdf[probdf$method=="neuralnet",5:8]))
print(ship_v_gigas)
pdf(file=bayesPlot_nn)
barplot(as.matrix(ship_v_gigas),beside=T,main=species)
dev.off()

ship_v_gigas = cbind(ship_pr=rowSums(probdf[probdf$method=="mnlogistic",1:4]),gigas_pr=rowSums(probdf[probdf$method=="mnlogistic",5:8]))
print(ship_v_gigas)
pdf(file=bayesPlot_mnlog)
barplot(as.matrix(ship_v_gigas),beside=T,main=species)
dev.off()

    
    probdf$species=species
    print("made it 2 in abc prior")
    write.table(file=probdffn,sep=",",row.names=F,rbind(probdf))

##
## build a random forest of the summary stats and train to predict the model indices
## run it on PCA transformed vars
##

if (!file.exists(rffn))
    {
        modelRF_pca = modelRF(pca.ref,cores=cores,modelColumn="introModel",group=list(ship=c("1","2"),gigas=c("3","4")))
        saveRDS(file=rffn,modelRF_pca)
    } else modelRF_pca=readRDS(file=rffn)

print (paste(date(),"ran pca RF"))

if (!file.exists(rfpfn))
    {
        modelRF_pca_pred = predRFmodel(modelRF_pca,pca.ref,cores=cores,modelColumn="introModel")
        saveRDS(file=rfpfn,modelRF_pca_pred)
    } else modelRF_pca_pred=readRDS(file=rfpfn)
print (paste(date(),"ran pca RFpred"))


if (!file.exists(rffnut))
    {
        modelRF = modelRF(untrans.ref,cores=cores,modelColumn="introModel",group=list(ship=c("1","2"),gigas=c("3","4")))
        saveRDS(file=rffnut,modelRF)
    } else modelRF=readRDS(file=rffnut)
print (paste(date(),"ran untrans RF"))


if (!file.exists(rfpfnut))
    {
        modelRF_pred = predRFmodel(modelRF,untrans.ref,cores=cores,modelColumn="introModel")
        saveRDS(file=rfpfnut,modelRF_pred)
    } else modelRF_pred=readRDS(file=rfpfnut)
    print (paste(date(),"ran untrans RFpred"))

    
if (diagnosticPlotRF)
{
pdf(file=paste0(abspath,"/figs/modelRF_pca_",priortag,".pdf"))
plotRFmodel(modelRF_pca,pca.ref)
dev.off()

pdf(file=paste0(abspath,"/figs/modelRF_pca_err_",priortag,".pdf"))
modelRF_pca_err = errRFmodel(modelRF_pca,pca.ref,cores=cores)
dev.off()

        pdf(file=paste0(abspath,"/figs/modelRF_",priortag,".pdf"))
        par(ask=FALSE)
        plotRFmodel(modelRF,untrans.ref)
        dev.off()
        
        pdf(file=paste0(abspath,"/figs/modelRF_err_",priortag,".pdf"))
        par(ask=FALSE)
        modelRF_err = errRFmodel(modelRF,untrans.ref,cores=cores)
        dev.off()


}

sink(rfsinkfn)
print(paste(species,date()))
print(modelRF_pca)
print(modelRF_pca_pred)
cat(paste0(species,"_RFPostProbPCA=",modelRF_pca_pred$post.prob,"\n"))
sink()

write.table(file=rfsumfn,sep=",",row.names=F,
            rbind(data.frame(species=species,type="pca",
                             chosen=modelRF_pca_pred$allocation,
                             post=modelRF_pca_pred$post.prob,oob=modelRF_pca$prior.err),
                  data.frame(species=species,type="untrans",
                             chosen=modelRF_pred$allocation,
                             post=modelRF_pred$post.prob,oob=modelRF$prior.err))
            )


    list(
        utfn=utfn,
        pcafn=pcafn,
        pcassfn=pcassfn,
        cvpfn=cvpfn,
        rffn=rffn,
        rfpfn=rfpfn,
        paramfn=paramfn,
        probdffn=probdffn,
        rfsumfn=rfsumfn,
        rfsinkfn=rfsinkfn,
        pcaplotfn=pcaplotfn,
        post_pca_fn=post_pca_fn,
        bayesPlot_nn=bayesPlot_nn,
        bayesPlot_mnlog=bayesPlot_mnlog
    )
}
                                        #============================================
#' Get Shipping Prior Probability
#'
#' This function reads a CSV file containing prior probabilities and filters
#' the probability associated with the current working directory.
#' It is intended to be used in a package eventually.
#'
#' @param fn A string specifying the file path to the CSV file containing prior probabilities.
#' @return Prints the working directory, the full priors dataframe, and the specific shipping probability.
#' @examples
#' getShipPrior("path/to/priors.csv")
#' @export
getShipPrior = function(fn)
{
    ##get the name of the executing directory (one up from src)
    dl = strsplit(getwd(),"/")[[1]]
    wd = dl[length(dl)-1]
    priors = read.csv(fn)
    probShip = priors$probShip[priors$simdir == wd]
    print(wd)
    print(priors)
    print(probShip)
}

#' Run ABC Group Replicate Analysis
#'
#' This function initializes the environment by sourcing data files and setting up
#' file paths. It then runs PCA transformations, builds and runs random forest models
#' on PCA-transformed and untransformed data, and outputs results.
#'
#' @param probShip Numeric, the prior probability of shipping (default is 0.5).
#' @param grouptag String, a tag for naming results from this group.
#' @param grouprep Integer, the replicate number.
#' @param datafn String, the filename of the R script that loads necessary data.
#' @param rfsample Integer, the number of samples for the random forest (default is 150000).
#' @param cores Integer, the number of cores to use for parallel processing (default is 1).
#' @param group List, additional grouping parameters for analysis.
#' @param verb Logical, whether to print verbose messages (default is TRUE).
#' @return A list of file names created and used during the process, as well as final model outcomes.
#' @examples
#' runABCgroupRep(probShip = 0.7, grouptag = "testGroup", grouprep = 1, datafn = "load_data.R")
runABCgroupRep = function(probShip = 0.5,     #prior prob of shipping 
                          grouptag = "grpcomp", #tag for naming results from this group
                          grouprep = 1,
                          datafn = "datafiles.R",
                          rfsample = 150000,
                          cores = 1,
                          group = list(),
                          verb = TRUE)
{
### initialize
    source(datafn)

    ##setup file names
    groupadd = paste0(grouptag,"_",grouprep)
    
    utfn=paste0(abspath,'/data/untransRef.RDS') #untrans ref file
    pcafn=paste0(abspath,'/data/pcaRef.RDS') #pca ref file
    rffn=paste0(abspath,'/data/modelRF_pca_',groupadd,'.RDS') #random forest on pca transformed data
    rfpfn=paste0(abspath,'/data/modelRF_pca_pred_',groupadd,'.RDS') ##pca model preds

    rffnut=paste0(abspath,'/data/modelRF_',groupadd,'.RDS')
    rfpfnut=paste0(abspath,'/data/modelRF_pred_',groupadd,'.RDS')


    rfsumfn=paste0(abspath,'/data/rfsum_',groupadd,'.csv') #write out posterir stuff for RF
    rfsinkfn=paste0(abspath,'/data/RFmodelOutput_',groupadd,'.txt') ##sink() outout for rf
    
    ##read in the raw summary stats, process and make a reference table
    genofile=paste(abspath,genofile,sep='/')
    mname=paste(abspath,mname,sep='/')



    
    ### start calcing!!


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
    print (paste(date(),"setup untrans ref"))

##########################
##
## pca transform the reference table and observations
##
if (!file.exists(pcafn))
    {
        pca.ref = make_pca_ref_table(untrans.ref,prop.variation=0.99)
        saveRDS(file=pcafn,pca.ref)
    } else pca.ref=readRDS(file=pcafn)
print (paste(date(),"setup pca ref"))

###
### make sure that the models adhere to the priors in the priorfile
###
### subset the models based on priors
###

print(paste0(date(),"about to run random forest"))

##
## build a random forest of the summary stats and train to predict the model indices
## run it on PCA transformed vars
##

modelRF_pca = modelRF(pca.ref,cores=cores,modelColumn="introModel",group=group)
print (paste(date(),"ran pca RF"))

    
modelRF_pca_pred = predRFmodel(modelRF_pca,pca.ref,cores=cores,modelColumn="introModel")
print (paste(date(),"ran pca RFpred"))


modelRF = modelRF(untrans.ref,cores=cores,modelColumn="introModel",group=group)
print (paste(date(),"ran untrans RF"))

modelRF_pred = predRFmodel(modelRF,untrans.ref,cores=cores,modelColumn="introModel")
print (paste(date(),"ran untrans RFpred"))


sink(rfsinkfn)
print(paste(species,date()))
print(modelRF_pca)
print(modelRF_pca_pred)
cat(paste0(species,"_RFPostProbPCA=",modelRF_pca_pred$post.prob,"\n"))
sink()


    rbind(data.frame(species=species,type="pca",
                     comparison=grouptag,
                     chosen=modelRF_pca_pred$allocation,
                     post=modelRF_pca_pred$post.prob,oob=modelRF_pca$prior.err),
          data.frame(species=species,type="untrans",
                     comparison=grouptag,
                     chosen=modelRF_pred$allocation,
                     post=modelRF_pred$post.prob,oob=modelRF$prior.err))
} ##end runABCgroupRep



