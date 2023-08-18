##
##
##
#' Plot Historic Events 
#'
#' This function visualizes migration events over time between different demes based on a simulation object.
#'
#' @param simObj a fsc object created with strataG tools
#' @param meta Metadata associated with the demes in the simulation.
#'
#' @details 
#' The function creates a plot with demes on the x-axis and time on the y-axis. Arrows indicate migration 
#' events between demes, with the direction of the arrow pointing from the source deme to the sink deme. 
#' The color of the arrow is determined by the source of the migration event, and the y-position indicates 
#' the time at which the migration event took place.
#'
#' The time axis is transformed using a logarithmic scale for better visualization.
#'
#' @return A plot visualizing the migration events over time.
#'
#' @examples 
#' \dontrun{
#' simObj <- some_simulation_object  # Replace with actual data or method to generate it
#' meta <- some_metadata_object      # Replace with actual data or method to generate it
#' plotevents(simObj, meta)
#' }
#'
#' @export
plotevents = function(simObj,meta)
{
    e = simObj$settings$events
    d = simObj$settings$demes
    d$index=0:(nrow(d)-1)
    meta=meta[order(meta$pop),]
    meta$deme=(1:nrow(meta))-1
    meta=meta[order(meta$source,meta$deme),]
    meta$idnum=1:nrow(meta)
    e=e[e$prop.migrants>0,]
    e=e[order(e$event.time),]
    
    e$event.time=log1p(e$event.time)
    lasts=rep(0,nrow(e))
    par(las=2)
    plot(1,type='n',xlim=range(c(e$source,e$sink)),ylim=c(0,max(e$event.time)),xlab="deme",ylab="time into past",axes=F)
    axis(1,at=sort(d$index),labels=paste(d$index,d$deme.name,sep="-"))
    axis(2)
    abline(v=sort(d$index),col="lightgray")
    last.tm=0
    for (tm in unique(e$event.time))
    {
        tms = e[e$event.time==tm,]
        sinks=unique(tms$sink)
        tms = tms[tms$prop.migrants>0,]
#        print(sinks)
        for (i in 1:nrow(tms))
        {
            
            lasts[tms$sink[i]+1] = tms$event.time[i]
            print(c(d$index[tms$source[i]],d$index[tms$sink[i]]))

            arrows(x0=tms$source[i],
                   x1=tms$sink[i],
                   y0=lasts[i],
                   y1=tms$event.time[i],
                   length=0.1,
                   col= as.numeric(unclass(factor(meta$source)[meta$deme==tms$source[i]])),
                   lwd=2
                   )
            
            
#            arrows(x0=meta$idnum[meta$deme==tms$source[i]],
#                   x1=meta$idnum[meta$deme==tms$sink[i]],
#                   y0=tms$event.time[i],
#                   y1=tms$event.time[i]
#                   )
#            abline(v=14)
        }
        for(s in sinks)
            arrows(x0=s,x1=s,y0=last.tm,y1=tm,length=0)
    last.tm=tm
    }
}
