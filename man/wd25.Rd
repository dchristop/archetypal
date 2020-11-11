\name{wd25}
\alias{wd25}
\docType{data}
\title{
2D data set created by 5 points for demonstration purposes
}
\description{
A data frame of 600 2D points 
}
\usage{data("wd25")}
\format{
matrix 600 x 2
}
\examples{
# Creation of data set "wd25" from 5 corner points:
set.seed(20191119)
p1 = c(3,2);p2 = c(4,6);p3 = c(7,8)
p4 = c(9,4);p5 = c(6,1)
dp = rbind(p1,p2,p3,p4,p5)
colnames(dp) = c('x','y')
pts=lapply(1:150, function(i,dp){
  c0 = runif(dim(dp)[1]);c0 = c0/sum(c0);pt0 = colSums(dp*c0)
  c1 = runif(3);c1 = c1/sum(c1);pt1 = colSums(dp[1:3,]*c1)
  c2 = runif(3);c2 = c2/sum(c2);pt2 = colSums(dp[c(4,5,1),]*c2)
  c3 = runif(3);c3 = c3/sum(c3);pt3 = colSums(dp[2:4,]*c3)
  rbind(pt0,pt1,pt2,pt3)
},dp)
df = do.call(rbind,pts)
rownames(df) = 1:dim(df)[1]
head(df)
# Check all equal
data("wd25")
all.equal(df,wd25)
# [1] TRUE
}
\keyword{datasets}
