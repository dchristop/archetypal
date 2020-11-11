\name{find_outmost_projected_convexhull_points}
\alias{find_outmost_projected_convexhull_points}
\title{Function which finds the outermost projected convex hull points
in order to be used as initial solution in archetypal analysis}
\usage{
find_outmost_projected_convexhull_points(df, kappas, npr = 2, rseed = NULL,
                                         doparallel = FALSE, nworkers = NULL,
                                         uniquerows = FALSE)
}
\arguments{
\item{df}{The n x d data frame that will be used for Archetypal Analysis}

\item{kappas}{The number of archetypes}

\item{npr}{The dimension of the projected subspaces. It can be npr = 1 (then there are d such subspaces), 
or npr > 1 (then we have C(d,npr) different subspaces)}

\item{rseed}{An integer to be used for the random seed if it will be necessary}

\item{doparallel}{If it is set to TRUE, then parallel processing will be performed.
That is absolutely required if n is very large and d>6.}

\item{nworkers}{The number of logical processors that will be used for computing the 
projected convex hulls, which they are always C(d,npr).}

\item{uniquerows}{If it is set to TRUE, then unique rows will be used for
computing distance matrix and less resources will be needed.}

}
\value{
A list with members:
\enumerate{
  \item outmost, the first kappas most frequent outermost points as rows of data frame
  \item outmostall, all the outermost points that have been found as rows of data frame
  \item outmostfrequency, a matrix with frequency and cumulative frequency for outermost rows
  \item usedrandom, an integer of randomly chosen rows, if it was necessary to complete the number of kappas rows
  \item chprojections, all the Convex Hulls of the different C(d,npr) projections, i.e. the coordinate projection subspaces 
  \item projected, a data frame with rows the unique points that have been projected in order to create 
  the relevant Convex Hulls of coordinate projection subspaces
}
}
\description{
Function which finds the outermost projected convex hull points
in order to be used as initial solution in archetypal analysis.
}
\details{
If \code{npr = 1}, then Convex Hull is identical with the range (\code{min},\code{max}) for the relevant variable, otherwise
the function uses the \code{\link[grDevices]{chull}} when \code{npr = 2} and the \code{\link[geometry]{convhulln}}
for \code{npr} > 2. See [1] and [2] respectively for more details.

First all available projections are being considered and their Convex Hull are being computed. Then either the unique (if \code{uniquerows = TRUE}) or all (if \code{uniquerows = FALSE}) associated data rows form a matrix and finally by using \code{\link[stats]{dist}} we find the kappas most frequent outermost rows. \cr
A special care is needed if the rows we have found are less than kappas. In that case, if a random sampling is necessary, 
the output \code{usedrandoms} informs us for the number of random rows and the \code{rseed} can be used for reproducibility.
}
\examples{
#
data("wd2") #2D demo 
df = wd2
yy = find_outmost_projected_convexhull_points(df, kappas = 3)
yy$outmost #the rows of 3 outmost projected convexhull points
yy$outmostall #all outmost found
yy$outmostfrequency #frequency table for all
yy$usedrandom #No random row was used
yy$chprojections #The Convex Hull of projection (one only here) 
yy$projected #the 9 unique points that created the one only CH
df[yy$outmost,] #the 3 outmost projected convexhull points
#
###
#
data("wd3") #3D demo 
df = wd3
yy = find_outmost_projected_convexhull_points(df, kappas = 4)
yy$outmost #the rows of 4 outmost projected convexhull points
yy$outmostall #all outmost found
yy$outmostfrequency #frequency table for all
yy$usedrandom #No random row was used
yy$chprojections #All the Convex Hulls of projections top coordinate planes
yy$projected #the 14 unique points that created all CHs
df[yy$outmost,] #the 4 outmost projected convexhull points
#
}
\references{
[1] Eddy, W. F. (1977). Algorithm 523: CONVEX, A new convex hull algorithm for planar sets.
ACM Transactions on Mathematical Software, 3, 411-412. doi: 10.1145/355759.355768.

[2] Barber, C.B., Dobkin, D.P., and Huhdanpraa, H.T., "The Quickhull algorithm for convex hulls"
ACM Trans. on Mathematical Software, 22(4):469-483, Dec 1996, http://www.qhull.org
}
\seealso{
\code{\link{find_furthestsum_points}}, \code{\link{find_outmost_convexhull_points}} 

\code{\link{find_outmost_partitioned_convexhull_points}} & \code{\link{find_outmost_points}}
}
