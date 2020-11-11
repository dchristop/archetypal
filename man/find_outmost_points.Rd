\name{find_outmost_points}
\alias{find_outmost_points}
\title{Function which finds the outermost points
in order to be used as initial solution in archetypal analysis}
\description{
Function which finds the outermost points
in order to be used as initial solution in archetypal analysis
}
\usage{
find_outmost_points(df, kappas)
}
\arguments{
\item{df}{The data frame with dimensions n x d}

\item{kappas}{The number of archetypes}
}
\value{
A list with members:
\enumerate{
  \item outmost, the first kappas most frequent outermost points as rows of data frame
  \item outmostall, all the outermost points that have been found as rows of data frame
  \item outmostfrequency, a matrix with frequency and cumulative frequency for outermost rows
}
}
\section{Warning}{
This is a rather naive way to find the outermost points of a data frame and 
it should be used with caution since for a n x d matrix we need in general
8 n^2/(2^30) GB RAM for numeric case. Check your machine and use it.
As a rule of thumb we advice its usage for n less or equal than 20000.
}
\examples{
data("wd2") #2D demo 
df = wd2
yy = find_outmost_points(df,kappas=3)
yy$outmost #the rows of 3 outmost points
yy$outmostall #all outmost found
yy$outmostfrequency #frequency table for all
df[yy$outmost,] #the 3 outmost points
#
###
#
data("wd3") #3D demo 
df = wd3
yy = find_outmost_points(df,kappas=4)
yy$outmost #the rows of 4 outmost points
yy$outmostall #all outmost found
yy$outmostfrequency #frequency table for all
df[yy$outmost,] #the 4 outmost points
}
\seealso{
\code{\link{find_furthestsum_points}},  \code{\link{find_outmost_convexhull_points}}, 

\code{\link{find_outmost_projected_convexhull_points}}, \cr

and \code{\link{find_outmost_partitioned_convexhull_points}}
}
