\name{find_outmost_partitioned_convexhull_points}
\alias{find_outmost_partitioned_convexhull_points}
\title{Function which finds the outmost convex hull points after making
np samples and findix convex hull for each of them.
To be used as initial solution in archetypal analysis.}
\description{
Function which finds the outmost convex hull points after making
np samples and findix convex hull for each of them.
To be used as initial solution in archetypal analysis
}
\usage{
find_outmost_partitioned_convexhull_points(df, kappas, np = 10,
  nworkers = NULL)
}
\arguments{
\item{df}{The data frame with dimensions n x d}

\item{kappas}{The number of archetypes}

\item{np}{The number of partitions that will be used (or the number of samples)}

\item{nworkers}{The number of logical processors that will be used}
}
\value{
A list with members:
\enumerate{
  \item outmost, the first kappas outmost points as rows of data frame
  \item outmostall, all the outmost points that have been found as rows of data frame
  \item outmostfrequency, a matrix with frequency and cumulative frequency for outmost rows
}
}
\examples{
data("wd2") #2D demo 
df=wd2
yy=find_outmost_partitioned_convexhull_points(df,kappas=3,nworkers = 1)
yy$outmost #the rows of 3 outmost points
df[yy$outmost,] #the 3 outmost points
yy$outmostall #all outmost rows
yy$outmostfrequency #their frequency
}
\seealso{
\code{\link{find_furthestsum_points}}, \code{\link{find_outmost_projected_convexhull_points}}, 

   \code{\link{find_outmost_convexhull_points}} & \code{\link{find_outmost_points}}
}
