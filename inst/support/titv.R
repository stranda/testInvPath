library(ape)

#Input: dat—an object of class ‘DNAbin’



titv<-function(dat){
  mat<-as.matrix(dat)
  res<-matrix(NA, ncol=dim(mat)[1], nrow=dim(mat)[1], dimnames=list(x=names(dat), y=names(dat)))
  for(i in 1:dim(mat)[1]){
    for(j in 1:dim(mat)[1])
      if (i>j) {
      vec<-as.numeric(mat[i,])+as.numeric(mat[j,])-8
      res[i,j]<-length(grep("200|56",vec)) #Transitions
      res[j,i]<-length(grep("152|168|88|104",vec)) #Transversions
    }
  }
  res
}


#Example

data(woodmouse)

ti<-titv(woodmouse)
tv<-t(ti)

tv[lower.tri(tv)] #Number of transversions
ti[lower.tri(ti)] #Number of transitions

#Saturation plot
dist<-dist.dna(woodmouse)

plot(ti[lower.tri(ti)]~dist)
points(tv[lower.tri(tv)]~dist, pch=20, col="red")
