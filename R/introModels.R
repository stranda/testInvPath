###
### different kinds of hard coded introduction models
###

#' shipping no longer than 30 days from sources
#' 
#' @param intros vector of names of intro regions
#' @param sources vector of names of source regions
#' @export
#' 
shipTbl30 <- function(intros=NULL,sources=NULL)
{
    
    m=structure(c(0L, 1L, 0L, 4L, 3L, 1L, 0L, 1L, 3L, 3L, 12L, 3L, 
                  3L, 6L, 11L, 6L, 6L, 2L, 17L, 23L),
                dim = 5:4,
                dimnames = list(
                    o = c("hok", "kag", "hon", "sea", "tok"),
                    d = c("EU", "NZ", "PNW", "Cali")), class = "table")
    if (is.null(intros)|is.null(sources))
        m
    else
    {
        tm=m[which(rownames(m)%in%sources),which(colnames(m)%in%intros)]
        if (length(intros)==1)
        {
            rn=names(tm)
            tm=matrix(tm,ncol=1)
            rownames(tm)=rn
            colnames(tm)=intros
            tm=as.table(tm)
        }
        tm
    }
}


#' shipping no longer than 60 days from sources
#' 
#' @param intros vector of names of intro regions
#' @param sources vector of names of source regions
#' @export
#' 

shipTbl60 <- function(intros=NULL,sources=NULL)
{
    m=structure(c(5L, 9L, 0L, 20L, 14L, 5L, 2L, 4L, 7L, 6L, 23L, 5L, 
                  3L, 11L, 22L, 10L, 11L, 2L, 25L, 35L),
                dim = 5:4,
                dimnames = list(
                    o = c("hok", "kag", "hon", "sea", "tok"),
                    d = c("EU",  "NZ", "PNW", "Cali")), class = "table")

    if (is.null(intros)|is.null(sources))
        m
    else     {
        tm=m[which(rownames(m)%in%sources),which(colnames(m)%in%intros)]
        if (length(intros)==1)
        {
            rn=names(tm)
            tm=matrix(tm,ncol=1)
            rownames(tm)=rn
            colnames(tm)=intros
            tm=as.table(tm)
        }
        tm
    }
}



#' gigas admixture-based introduction model
#' 
#' @param intros vector of names of intro regions
#' @param sources vector of names of source regions
#' @export
#' 

gigasTblAdmix = function(intros=NULL,sources=NULL)
{
    m=matrix(
###   Ar  Ch  EN    ES    NZ    PNW   Cali
    c(0.5,0.5,0.5,  0.1,   0,    0.5,    0,  #hok
      0.5,0.5,0.5,  0.9,   0,    0.5, 0.20,  #hon
      0,  0,  0,    0,   0.5,    0,   0.40,  #tok
      0,  0,  0,    0,   0.5,    0,   0.40,  #sea
      0,  0,  0,    0,     0,    0,      0,  #kag
      0,  0,  0,    0,     0,    0,      0), #nonSource
    byrow=T,ncol=7)

    colnames(m)=c("Ar","Ch","EN","ES","NZ","PNW","Cali")
    rownames(m)=c("hok","hon","tok","sea","kag","nonSource")

    m=as.table(cbind(m,EU=rowSums(m[,3:4])))
    if (is.null(intros)|is.null(sources))
        m
    else     {
        tm=m[which(rownames(m)%in%sources),which(colnames(m)%in%intros)]
        if (length(intros)==1)
        {
            rn=names(tm)
            tm=matrix(tm,ncol=1)
            rownames(tm)=rn
            colnames(tm)=intros
            tm=as.table(tm)
        }
        tm
    }
}



#' gigas machine learning model based introduction table 
#' 
#' @param intros vector of names of intro regions
#' @param sources vector of names of source regions
#' @export
#' 

gigasTblML <-function(intros=NULL,sources=NULL)
{
    m = structure(c(0L, 34L, 0L, 0L, 0L, 0L, 1L, 14L, 0L, 3L, 0L, 0L, 
1L, 68L, 0L, 0L, 0L, 0L, 3L, 58L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 
55L, 0L, 0L, 2L, 59L, 0L, 0L, 0L, 0L, 0L, 36L, 0L, 37L, 0L, 2L
), dim = 6:7, dimnames = list(pred = c("hok", "hon", "tok", "sea", 
"kag", "nonSource"), actual = c("Arg", "Chile", "EuropeNorth", 
"EuropeSouth", "NZ", "PNW", "Cali")), class = "table")

m=cbind(m,EU=rowSums(m[,3:4]))
if (is.null(intros)|is.null(sources))
    m
else     {
        tm=m[which(rownames(m)%in%sources),which(colnames(m)%in%intros)]
        if (length(intros)==1)
        {
            rn=names(tm)
            tm=matrix(tm,ncol=1)
            rownames(tm)=rn
            colnames(tm)=intros
            tm=as.table(tm)
        }
        tm
    }
}

#' shipping no longer than 30 days from sources
#' 
#' @param intros vector of names of intro regions
#' @param sources vector of names of source regions
#' @export
#' 

equalSourceTbl <- function(intros=NULL,sources=NULL)
{
    
    m=structure(rep(1,20),
                dim = 5:4,
                dimnames = list(
                    o = c("hok", "kag", "hon", "sea", "tok"),
                    d = c("EU", "NZ", "PNW", "Cali")), class = "table")
    if (is.null(intros)|is.null(sources))
        m
    else     {
        tm=m[which(rownames(m)%in%sources),which(colnames(m)%in%intros)]
        if (length(intros)==1)
        {
            rn=names(tm)
            tm=matrix(tm,ncol=1)
            rownames(tm)=rn
            colnames(tm)=intros
            tm=as.table(tm)
        }
        tm
    }
}

