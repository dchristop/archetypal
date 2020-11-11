\name{find_furthestsum_points}
\alias{find_furthestsum_points}
\title{Function which finds the furthest sum points
in order to be used as initial solution in archetypal analysis}
\description{
Function which finds the furthest sum points
in order to be used as initial solution in archetypal analysis.
}
\usage{
find_furthestsum_points(df, kappas, nfurthest = 100, nworkers = NULL,
                        sortrows = TRUE, doparallel = TRUE)
}
\arguments{
\item{df}{The data frame with dimensions n x d}

\item{kappas}{The number of archetypes}

\item{nfurthest}{The number of applications for FurthestSum algorithm}

\item{nworkers}{The number of logical processors that will be used.
Hint: set it such that \code{nfurthest} can be an exact multiple of \code{nworkers}.}

\item{sortrows}{If it is TRUE, then rows will be sorted}

\item{doparallel}{If it is set to TRUE, then parallel processing will be performed
for the \code{nfurthest} applications of algorithm}
}
\value{
A list with members:
\enumerate{
  \item outmost, the first kappas furthest sum points as rows of data frame
  \item outmostall, all the furthest sum points that have been found as rows of data frame
  \item outmostfrequency, a matrix with frequency and cumulative frequency for furthest sum rows
}
}
\examples{
data("wd3") #3D demo 
df = wd3
yy = find_furthestsum_points(df, kappas = 4, nfurthest = 10, nworkers = 2)
yy$outmost
yy$outmostall
yy$outmostfrequency
}
\seealso{
\code{\link{FurthestSum}}
}
