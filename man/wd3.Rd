\name{wd3}
\alias{wd3}
\docType{data}
\title{
3D data set for demonstration purposes
}
\description{
A data frame of 100 3D points 
}
\usage{data("wd3")}
\format{
matrix 100 x 3
}
\examples{
# Creation of data set "wd3" from 4 corner points:
p1 = c(3,0,0);p2 = c(0,5,0)
p3 = c(3,5,7);p4 = c(0,0,0)
# The data frame of generators
dp = data.frame(rbind(p1,p2,p3,p4))
colnames(dp) = c("x","y","z")
dp = dp[chull(dp),]
set.seed(9102) 
df = data.frame(t(sapply(1:100, function(i,dp){
  cc = runif(4)
  cc = cc/sum(cc) 
  colSums(dp*cc)
},dp)))
colnames(df) = c("x","y","z")
head(df)
# Check all.equal to "wd3"
data(wd3)
all.equal(df,wd3)
# [1] TRUE
}
\keyword{datasets}
