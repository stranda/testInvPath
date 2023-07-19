#' See if a gtypes and a meta object have concordant information
#' @param meta  meta data frame
#' @param gin gtypes object (must have populations as strata).
#'
#' @description Tests to see if the populations in meta match up to the strata in gin
#' 
#' @return boolean T if all required cols are present in meta
#'
#' @export
concordantGtypesMeta = function(meta,gin)
{
    ok=T

    if (!is.meta(meta))
    {
        print("wrong columns in meta")
        ok = F
    }

    if (class(gin)!="gtypes")
    {
        print("gin is not a strataG gtypes object")
        ok=F
    }
    
    if (ok) #make sure that there are the same strata and pops
    {
        s=unique(getStrata(gin))
        p=unique(meta$pop)
        if (sum(s%in%p)!=length(s))
        {
            print(paste("these strata are missing from meta",paste(s[!s%in%p], collapse=", ")))
            ok=F
        }
        if (sum(p%in%s)!=length(p))
        {
            print(paste("these pops in meta are missing from gin strata",paste(p[!p%in%s], collapse=", ")))
            ok=F
        }
    }
       
    ok
}



#' Test if meta has the correct columns
#' @param meta meta data
#' @param required a vector of names of columns required in meta
#'
#' @returns boolean T if all required cols are present in meta
#' 
#' @export
is.meta = function(meta,required=c("pop","longpop","source","reg"))
{
    if (is.null(required))
        required=c("pop","longpop","source","reg")

    if (sum(required%in%names(meta))==length(required))
        TRUE
    else
        FALSE
}
