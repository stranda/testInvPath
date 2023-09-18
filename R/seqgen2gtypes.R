#' Convert take a tree create sequences using phyclust::seqgen and convert to gtypes
#'
#' This function generates sequences based on a given phylogenetic tree and parameters 
#' using the seqgen function from the phyclust package. It then converts these sequences 
#' into gtype format using the strataG::sequence2gtypes function.
#'
#' @param tree A phylogenetic tree to base sequence generation on.
#' @param params A list containing parameters for sequence generation.
#' @param temp.file A temporary file location to write sequences to.
#' @param demes A list or data frame containing information about demes or populations.
#' @param mu A numeric parameter for mutation rate. Default is 1.0.
#'
#' @return A gtypes object.
#' @export
#'
seqgen2gtype = function(tree, params , temp.file, demes,mu=1.0)
{
    opts=paste0("-q -mHKY -d",mu," ",
                "-l",ncol(params$gin@sequences@dna$gene.1)#seq length
                )
    print(opts)
    phyclust::seqgen(opts,tree,temp.file=temp.file)
    db=read.dna(temp.file,format="sequential")

    rownames(db)=paste(sapply(strsplit(rownames(db),"\\."),function(x)as.numeric(x[2])),
          sapply(strsplit(rownames(db),"\\."),function(x)as.numeric(x[1])),
          sep="_")
    
    a=sequence2gtypes(db)
    pops=demes$deme.name[sapply(strsplit(getIndNames(a),"_"),function(x){as.numeric(x[1])})]

    names(pops)=getIndNames(a)
    setStrata(a)=pops
    a #return the gtpyes object
}

