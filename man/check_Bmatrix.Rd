\name{check_Bmatrix}
\alias{check_Bmatrix}
\title{Function which checks B matrix of Archetypal Analysis  Y ~ A B Y
in order to find the used rows for creating each archetype
and the relevant used weights.}
\description{
Function which checks B matrix of Archetypal Analysis  Y ~ A B Y
in order to find the used rows for creating each archetype
and the relevant used weights.
}
\usage{
check_Bmatrix(B, print.details = TRUE)
}
\arguments{
\item{B}{The \eqn{kappas \times n} matrix such that Y ~ Y - ABY or Frobenius norm ||Y-ABY|| is minimum}

\item{print.details}{If set to TRUE, then results are printed out.}
}
\value{
A list with members:
\enumerate{
  \item used_rows, a list with used rows for creating each archetype
  \item used_weights, a list with the relevant weights that have been used
}
}
\examples{
# Create a small 2D data set from 3 corner-points:
p1=c(1,2);p2=c(3,5);p3=c(7,3) 
dp=rbind(p1,p2,p3);dp
set.seed(916070)
pts=t(sapply(1:20, function(i,dp){
  cc=runif(3)
  cc=cc/sum(cc)
  colSums(dp*cc)
},dp))
df=data.frame(pts)
colnames(df)=c("x","y")
# Run AA:
aa=archetypal(df=df,kappas = 3,verbose=FALSE)
# Check B matrix:
B=aa$B
yy=check_Bmatrix(B)
yy$used_rows
yy$used_weights
#
}
\seealso{
\code{\link{archetypal}}
}
