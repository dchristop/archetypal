\name{wd2}
\alias{wd2}
\docType{data}
\title{
2D data set for demonstration purposes
}
\description{
A data frame of 100 2D points
}
\usage{data("wd2")}
\format{
matrix 100 x 2
}
\examples{
# Creation of data set "wd2" from 3 corner-points:
p1 = c(1,2);p2 = c(3,5);p3 = c(7,3) 
dp = rbind(p1,p2,p3);dp
set.seed(9102)
pts = t(sapply(1:100, function(i,dp){
  cc = runif(3)
  cc = cc/sum(cc)
  colSums(dp*cc)
},dp))
df = data.frame(pts)
colnames(df) = c("x","y")
head(df)
# Check all equal:
data(wd2)
all.equal(wd2,df)
# [1] TRUE
}
\keyword{datasets}
