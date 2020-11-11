\name{find_outmost_partitioned_convexhull_points}
\alias{find_outmost_partitioned_convexhull_points}
\title{Function which finds the outermost convex hull points after making
np samples and finding convex hull for each of them.
}
\description{
Function which finds the outermost convex hull points after making
\code{np} samples and finding convex hull for each of them.
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
  \item outmost, the first kappas most frequent outermost points as rows of data frame
  \item outmostall, all the outermost points that have been found as rows of data frame
  \item outmostfrequency, a matrix with frequency and cumulative frequency for outermost rows
}
}
\examples{
data("wd2") #2D demo 
df = wd2
yy = find_outmost_partitioned_convexhull_points(df, kappas = 3, nworkers = 2)
yy$outmost #the rows of 3 outermost points
df[yy$outmost,] #the 3 outermost points
yy$outmostall #all outermost rows
yy$outmostfrequency #their frequency
}
\seealso{
\code{\link{find_furthestsum_points}}, \code{\link{find_outmost_projected_convexhull_points}}, 

   \code{\link{find_outmost_convexhull_points}} & \code{\link{find_outmost_points}}
}
