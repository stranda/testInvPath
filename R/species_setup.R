#' Set up the file hierarchy and data files
#'
#' @param root  name of the root directory in the hierarchy to get information
#' @param fas name of a fasta file containing haplotypes. will be converted to a gtypes
#' @param genofile filename of a gtypes object in .RDS format (must have populations as strata). Either fas or genofilehave to be specified
#' @param mname name of file containing metadata that links pops to strata
#' @param newdir path to the new directory to set things up for a species
#' @param species the name of the species under consideration
#' @param dataType type of data ("sequence", "snp", "microsatellite")
#' @param ploidy if microsatellite always set to 2.  Otherwise 1 unless overridden.  for sequences, the adjustment is ploidy*sample.size for the number sampled per deme.
#' @param nativeTopology a named matrix describing the topology and relative times for the history of a species in the native range.  If null defaults to output of function nativeHistory()
#' @param fsc_exec executable for fastsimcoal
#' @export
species_setup <- function(root="test",fas=NULL,indmeta=NULL,genofile=NULL,mname=NULL,newdir=NULL,dataType="sequence",species="test",nativeTopology=NULL,fsc_exec="fsc27",popPairwise=FALSE,use.seqgen=FALSE,ploidy=1)
{
    if (is.null(fas)&is.null(genofile)) stop("there has to be some genetic data")
    if (is.null(mname)) stop("the genetic data must be accompanied by metadata linking populatins to strata")

    if (dataType=="microsatellite") ploidy=2

    if (!is.null(fas))
    {

        if (!file.exists(fas)){
            if (file.exists(paste0(root,"/",fas)))
                fas=paste0(root,"/",fas)
            else if (file.exists(paste0(root,"/data/",fas)))
                fas=paste0(root,"/data/",fas)
            else stop(paste("cant find fasta file ",fas))
        }
        if (!file.exists(indmeta)){
                        if (file.exists(paste0(root,"/",indmeta)))
                            indmeta=paste0(root,"/",indmeta)
                        else if (file.exists(paste0(root,"/data/",indmeta)))
                            indmeta=paste0(root,"/data/",indmeta)
                        else stop(paste("cant find fasta file ",indmeta))                    
        }
#        print("running fasta2gin")
        g=fasta2gin(fas,indmeta)
        genofile=paste0(root,"/","gtype_from_fasta.RDS")
        saveRDS(file=paste0(genofile),g)
    }
    
    if (!is.null(genofile))
        if (!file.exists(genofile)){
            if (file.exists(paste0(root,"/",genofile)))
                genofile=paste0(root,"/",genofile)
            else   if (file.exists(paste0(root,"/data/",genofile)))
                genofile=paste0(root,"/data/",genofile)
            else stop(paste("cant find genofile file ",genofile))
        }

    if (!file.exists(mname))
    {
        if (file.exists(paste0(root,"/",mname)))
            mname=paste0(root,"/",mname)
        else if (file.exists(paste0(root,"/data/",mname)))
            mname=paste0(mname)
        else stop(paste("cant find meta file ",mname))
    }

    meta1 = read.csv(mname)
    if (all(c("pop", "source") %in% names(meta1)))
    {
        if (!any(meta1$source %in% colnames(gigasTblML())))
        {
            print(paste("metafile sources: ",paste(unique(meta1$source),collapse=", ")))
            print(paste("colnames: ",paste(unique(colnames(gigasTblML())),collapse=", ")))
            stop("the sources in the metafile are not found in the colnames of the the introduction matrices")}
        
    } else stop("problems with metafile, missing columns pop and source")

    gt=gigasTblML()
    itrs = unique(colnames(gt)[colnames(gt)%in%meta1$source])
    srcs = unique(rownames(gt)[rownames(gt)%in%meta1$source])

    if (!dir.exists(newdir))
    {
        dir.create(newdir)
        if (!dir.exists(paste0(newdir,"/data"))) dir.create(paste0(newdir,"/data"))
        if (!dir.exists(paste0(newdir,"/src"))) dir.create(paste0(newdir,"/src"))
        if (!dir.exists(paste0(newdir,"/results"))) dir.create(paste0(newdir,"/results"))
        if (!dir.exists(paste0(newdir,"/logs"))) dir.create(paste0(newdir,"/logs"))
        if (!dir.exists(paste0(newdir,"/gathered"))) dir.create(paste0(newdir,"/gathered"))
        if (!dir.exists(paste0(newdir,"/archivedCSV"))) dir.create(paste0(newdir,"/archivedCSV"))
        if (!dir.exists(paste0(newdir,"/figs"))) dir.create(paste0(newdir,"/figs"))
    
    file.copy(mname,paste0(newdir,"/data/metadata.csv"))
    file.copy(mname,paste0(newdir,"/data/"))
    file.copy(genofile,paste0(newdir,"/data/gin.RDS"))
    file.copy(genofile,paste0(newdir,"/data/"))

    
    if (!is.null(fas)) {
        file.copy(fas,paste0(newdir,"/data/"))
        file.copy(indmeta,paste0(newdir,"/data/"))
    }

    
############################  copy some R scripts 
    file.copy(paste0(system.file(package="testInvPath","skeletons"),"/runReps.R"),paste0(newdir,"/src/"))
        file.copy(paste0(system.file(package="testInvPath","skeletons"),"/abc.R"),paste0(newdir,"/src/"))
############################ copy shell Scripts
    file.copy(paste0(system.file(package="testInvPath","skeletons"),"/runSimulation.sh"),paste0(newdir,"/"))
    Sys.chmod(paste0(newdir,"/runSimulation.sh"),mode="755")
    file.copy(paste0(system.file(package="testInvPath","skeletons"),"/gatherCSV.sh"),paste0(newdir,"/"))
    Sys.chmod(paste0(newdir,"/gatherCSV.sh"),mode="755")
    file.copy(paste0(system.file(package="testInvPath","skeletons"),"/runABC.sh"),paste0(newdir,"/"))
    Sys.chmod(paste0(newdir,"/runABC.sh"),mode="755")
################### create some files

    ########datafiles file
    cat(file=paste0(newdir,"/src/datafiles.R"),append=F,"#These are the relative paths to the data files\n\n")
    if (!is.null(fas)){

        cat(file=paste0(newdir,"/src/datafiles.R"),append=T,paste0("fas='data/",basename(fas),"' #fasta source\n\n"))
        cat(file=paste0(newdir,"/src/datafiles.R"),append=T,paste0("indmeta='data/",basename(indmeta),"' #individual-level metadata (strata assignment, if needed)\n\n"))

    } else {
        cat(file=paste0(newdir,"/src/datafiles.R"),append=T,paste0("fas=NULL  #fasta source (NULL=none)\n"))
        cat(file=paste0(newdir,"/src/datafiles.R"),append=T,paste0("indmeta=NULL\n\n")) 

        }

    cat(file=paste0(newdir,"/src/datafiles.R"),append=T,paste0("genofile='data/",basename(genofile),"' #RDS file with gtypes inside\n\n"))
    cat(file=paste0(newdir,"/src/datafiles.R"),append=T,paste0("mname='data/",basename(mname),"' #name of the metafile for  pops\n\n"))

    cat(file=paste0(newdir,"/src/datafiles.R"),append=T,paste0("abspath='",newdir,"' # absolute path to species, relative paths above are within this parent directory.  Change this if needed\n\n"))


    cat(file=paste0(newdir,"/src/datafiles.R"),append=T,paste0("######## sources and introductions\n\n"))

    cat(file=paste0(newdir,"/src/datafiles.R"),
        append=T,
        paste0("sources=",capture.output(dput(srcs,file=""))," # sources in the metadata and in the intro matrices\n\n"))

        cat(file=paste0(newdir,"/src/datafiles.R"),append=T,paste0("sources=sources[sources!='nonSource'] # comment out to include nonSource\n\n"))


    cat(file=paste0(newdir,"/src/datafiles.R"),append=T,paste0("intros=",capture.output(dput(itrs,file=""))," # introductions in the metadata and in the intro matrices\n\n"))

    if (is.null(nativeTopology)) nativeTopology=nativeHistory(sources=srcs)     

    cat(file=paste0(newdir,"/src/datafiles.R"),append=T,paste0("# native range topology:\n"))

    cat(file=paste0(newdir,"/src/datafiles.R"),append=T,"nativeTopology=",capture.output(dput(nativeTopology,file="")),"\n\n")
    
        cat(file=paste0(newdir,"/src/datafiles.R"),append=T,paste0("dataType='",dataType,"' #type of genetic data\n\n"))

        cat(file=paste0(newdir,"/src/datafiles.R"),append=T,paste0("ploidy='",ploidy,"' #ploidy of loci\n\n"))

        cat(file=paste0(newdir,"/src/datafiles.R"),append=T,paste0("popPairwise=",popPairwise," #caculate pairwise _population_ statistics (regional are retained)\n\n"))

        cat(file=paste0(newdir,"/src/datafiles.R"),append=T,paste0("use.seqgen=",use.seqgen," #fsc makes a tree but use seqgen to add sequences\n\n"))
        
    cat(file=paste0(newdir,"/src/datafiles.R"),append=T,paste0("species='",species,"' #species name\n\n"))

        cat(file=paste0(newdir,"/src/datafiles.R"),append=T,paste0("fsc_exec='",fsc_exec,"' #fastsimcoal executable name\n\n"))
    } else { #Species directory structure already exists
        print(paste0(newdir," already exists.  To be safe I didn't overwrite anything.  Delete or rename it"))
    }
}    
    
    
