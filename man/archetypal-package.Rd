\docType{package}
\name{archetypal-package}
\alias{archetypal-package}
\title{Finds the Archetypal Analysis of a Data Frame}
\description{
Performs archetypal analysis by using Convex Hull approximation under a full control of all algorithmic parameters. 
It contains functions useful for finding the proper initial approximation, the optimal number of archetypes and 
function for applying main algorithm.
}
\section{Compute Archerypal Analysis (AA)}{
The main function is \code{\link{archetypal}} which is a modification of PCHA algorithm, see [1], [2],
suitable for R language. It provides control to the entrire set of involved parameters  and has two main options:   
\enumerate{
 \item initialrows=NULL, then a method from 'projected_convexhull', 'convexhull', 'partitioned_convexhull',
'furthestsum', 'outmost', 'random' is used
 \item initialrows={a vector of kappas rows}, then given rows form the initial solution for AA
} 
This is the main function of the package, but extensive trials has shown that:
\itemize{
 \item AA may be very difficult to run if a random initial solution has beeen chosen
 \item for the same data set the final Sum of Squared Errors (SSE) may be much smaller
 if initial solution is close to the final one
 \item even the quality of AA done is affected from the starting point
 }
This is the reason why we have developed a whole set of methods for choosing initial solutionfor main PCHA algorithm.
} 
\section{Find a time efficient initial approximation for AA}{ 
There are three functions that work with the Convex Hull (CH) of data set.
 \enumerate{
    \item  \code{\link{find_outmost_convexhull_points}} computes the CH of all points    
    \item \code{\link{find_outmost_projected_convexhull_points}} computes the CH for
    all possible combinations of variables taken by n (default=2)
    \item \code{\link{find_outmost_partitioned_convexhull_points}} makes partitions
   of data frame, then computes CH for each partition and finally gives the CH of overall union  
  }
   
 The most simple method for estimating an initial solution is \code{\link{find_outmost_points}}
 where we just compute the "outmost points", ie those that are the most frequent outmost for all
 available points.
   
 The default method "FurthestSum" (FS) of PCHA (see [1], [2]) is used by \code{\link{find_furthestsum_points}} which applies
 FS for many times (default=10) and then finds the most frequent points. 
 
 Of course "random" method is available for comparison reasons and that gives a random
 set of kappas points as initial solution. 
 
 All methods give the number of rows for the input data frame.
}
\section{Find the optimal number of archetypes (kappas)}{
For that task \code{\link{find_optimal_kappas}} is available which 
runs for each kappas from 1 to maxKappas (default=15) ntrials (desault=10) times AA, 
stores SSE, VarianceExplained from each run and then computes knee point by nusing UIK method, see [3].
} 
\section{Evaluate the quality of Archerypal Analysis (AA)}{
By using function \code{\link{check_Bmatrix}} we can evaluate the overall quality of 
applied method and algorithm. Quality can be considered high:
\enumerate{
  \item if every archetype is being created by a small number of points
  \item if relevant weights are not numerically insignificant
}
Of course we must take into account the SSE and VarianceExplained, but if we have to compare two solutions
with similar termination status, then we must choose that of the simplest B matrix form.
} 
\note{
Bug reports and feature requests can be sent to
  \email{dchristop@econ.uoa.gr}.
}
\author{
\strong{Maintainer}: Demetris Christopoulos \email{dchristop@econ.uoa.gr}

Other contributors:
\itemize{
  \item David Midgley \email{david.midgley@insead.edu} [contributor]
  \item INSEAD Fontainebleau France [funder, copyright holder]
}
}
\references{
{
[1] M Morup and LK Hansen, "Archetypal analysis for machine learning and data mining", Neurocomputing (Elsevier, 2012).
 https://doi.org/10.1016/j.neucom.2011.06.033.

[2] Source: http://www.mortenmorup.dk/index_files/Page327.htm , last accessed 2019-06-07

[3] Christopoulos, Demetris T., Introducing Unit Invariant Knee (UIK) As an Objective Choice for
 Elbow Point in Multivariate Data Analysis Techniques (March 1, 2016). 
Available at SSRN: https://ssrn.com/abstract=3043076 or http://dx.doi.org/10.2139/ssrn.3043076

}
}
\seealso{
\code{\link{archetypal}}
}
\keyword{archetypal}
\keyword{PCHA}
\keyword{convex hull}