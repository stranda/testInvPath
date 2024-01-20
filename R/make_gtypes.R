###
### make a gtypes object from fasta data and some meta data
###

fasta2gin = function(fas=NULL,indmeta=NULL)
    {
###        fas = "../data/Uwai_141inds_2loc.fas"
###        mn =  "../data/Uwai_141inds.meta.csv"
        mn=indmeta
        lmeta = read.csv(mn)
#        print(head(lmeta))
        rownames(lmeta)=lmeta$indID
        fasdat=strataG::read.fasta(fas)
#print("converting fasta to gtype in the next line")        
        g = strataG::sequence2gtypes(fasdat)
        
#        print(getIndNames(g))

        sv = lmeta$popID
        names(sv)=rownames(lmeta)
        strataG::setStrata(g) = sv

        g

    }
