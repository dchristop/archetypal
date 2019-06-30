\name{find_outmost_convexhull_points}
\alias{find_outmost_convexhull_points}
\title{Function which finds the outmost convex hull points
in order to be used as initial solution in archetypal analysis}
\usage{
find_outmost_convexhull_points(df, kappas)
}
\arguments{
\item{df}{The data frame with dimensions n x d}
\item{kappas}{The number of archetypes}
}
\value{
A list with members:
\enumerate{
  \item outmost, the first kappas outmost points as rows of data frame
  \item outmostall, all the outmost points that have been found as rows of data frame
  \item outmostfrequency, a matrix with frequency and cumulative frequency for outmost rows
}
}
\description{
Function which finds the outmost convex hull points
in order to be used as initial solution in archetypal analysis
}
\details{
This function uses the \code{\link[grDevices]{chull}} when n=2 (see [1], [2]) and the \code{\link[geometry]{convhulln}}
for n>2 (see [3])  cases.
}
\examples{
data("wd2") #2D demo 
df=wd2
yy=find_outmost_convexhull_points(df,kappas=3)
yy$outmost #the rows of 3 outmost points
df[yy$outmost,] #the 3 outmost points
yy$outmostall #all outmost cH rows
yy$outmostfrequency #their frequency
#
###
#
data("wd3") #3D demo 
df=wd3
yy=find_outmost_convexhull_points(df,kappas=4)
yy$outmost #the rows of 4 outmost points
df[yy$outmost,] #the 4 outmost points
yy$outmostall #all outmost cH rows
yy$outmostfrequency #their frequency
}
\references{
[1] Eddy, W. F. (1977). A new convex hull algorithm for planar sets.
  ACM Transactions on Mathematical Software, 3, 398-403. doi: 10.1145/355759.355766.

[2] Eddy, W. F. (1977). Algorithm 523: CONVEX, A new convex hull algorithm for planar sets [Z].
  ACM Transactions on Mathematical Software, 3, 411-412. doi: 10.1145/355759.355768.

[3] Barber, C.B., Dobkin, D.P., and Huhdanpaa, H.T., "The Quickhull algorithm for convex hulls"
  ACM Trans. on Mathematical Software, 22(4):469-483, Dec 1996, http://www.qhull.org
}
\seealso{
\code{\link{find_furthestsum_points}}, \code{\link{find_outmost_projected_convexhull_points}},  

   \code{\link{find_outmost_partitioned_convexhull_points}} &  \code{\link{find_outmost_points}}
}
