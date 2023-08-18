#' set up the history for the native range
#' @param sources vector of source regions (should be same as rownames for intro tables)
#'
#'the native range topology maps the times t1-t3 to the 4 regions in japan
#'will use this in the generate coalescence model functions to determine within
#'native range coalescence.  This works by going down columns of the matrix.  Hokkaido
#'coalesces into honshu at time 1 and honshu into tokyo at time 3 (at least for the
#'original gvermSNP file.  if this is downstream, regions likely differ
#' @export
nativeHistory = function(sources=c("kag","hon","sea","tok","hok"),
                                            #k,h,s,t,hk   
                         topology= matrix(c(0,0,0,0,0,  #kag
                                            0,0,0,0,1,  #hon
                                            1,0,0,0,0,  #sea
                                            0,3,2,0,0,  #tok
                                            0,0,0,0,0), #hok
                                          byrow=T,nrow=5,
                                          dimnames=list(c("kag","hon","sea","tok","hok"),
                                                        c("kag","hon","sea","tok","hok"))
                                          )
                               )
{
    newtop=topology[which(rownames(topology)%in%sources),which(colnames(topology)%in%sources)]
    if (min(colSums(newtop)>0)) {stop("no root for native topology, need a column in the matrix with all zeros")}
    if (sum(colSums(newtop)==0)>1) {
        stop(paste("there are multiple places in the native topology where lineages do not connect:",paste(colnames(newtop)[colSums(newtop)==0],collapse=" ,")))}
    newtop
}
