\docType{package}
\name{archetypal-package}
\alias{archetypal-package}
\title{Finds the Archetypal Analysis of a Data Frame}
\description{
Performs archetypal analysis by using Principal Convex Hull Analysis (PCHA) 
under a full control of all algorithmic parameters. 
It contains a set of functions for determining the initial solution, the optimal 
algorithmic parameters and the optimal number of archetypes.
Post run tools are also available for the assessment of the derived solution. 
}
\section{Compute Archetypal Analysis (AA)}{
The main function is \code{\link{archetypal}} which is a variant of PCHA algorithm, see [1], [2],
suitable for R language. It provides control to the entire set of involved parameters  and has two main options:   
\enumerate{
 \item initialrows = NULL, then a method from "projected_convexhull", "convexhull",\cr
 "partitioned_convexhul", "furthestsum", "outmost", "random" is used
 \item initialrows = (a vector of kappas rows), then given rows form the initial solution for AA
} 
This is the main function of the package, but extensive trials has shown that:
\itemize{
 \item AA may be very difficult to run if a random initial solution has been chosen
 \item for the same data set the final Sum of Squared Errors (SSE) may be much smaller
 if initial solution is close to the final one
 \item even the quality of AA done is affected from the starting point
 }
This is the reason why we have developed a whole set of methods for choosing initial solution for the PCHA algorithm.
} 
\section{Find a time efficient initial approximation for AA}{ 
There are three functions that work with the Convex Hull (CH) of data set.
 \enumerate{
    \item  \code{\link{find_outmost_convexhull_points}} computes the CH of all points    
    \item \code{\link{find_outmost_projected_convexhull_points}} computes the CH for
    all possible combinations of variables taken by \code{npr} (default=2)
    \item \code{\link{find_outmost_partitioned_convexhull_points}} makes \code{np} partitions
   of data frame (defualt=10), then computes CH for each partition and finally gives the CH of overall union  
  }
   
 The most simple method for estimating an initial solution is \code{\link{find_outmost_points}}
 where we just compute the outermost points, i.e. those that are the most frequent outermost for all
 available points.
   
 The default method "FurthestSum" (FS) of PCHA (see [1], [2]) is used by \code{\link{find_furthestsum_points}} which applies
 FS for \code{nfurthest} times (default=10) and then finds the most frequent points. 
 
 Of course "random" method is available for comparison reasons and that gives a random
 set of kappas points as initial solution. 
 
 All methods give the number of rows for the input data frame as integers. Attention needed if your data frame
 has row names which are integers but not identical to \code{1:dim(df)[1]}.
}
\section{Find the optimal number of archetypes}{
For that task \code{\link{find_optimal_kappas}} is available which 
runs for each kappas from 1 to maxkappas (default=15) ntrials (default=10) times AA, 
stores SSE, VarianceExplained from each run and then computes knee or elbow point by using UIK method, see [3].
} 
\section{Determining the optimal updating parameters}{
Extensive trials have shown us that choosing the proper values for algorithmic updating parameters 
(\code{muAup, muAdown, muBup, muBdown}) can speed up remarkably the process. That is the task of
 \code{\link{find_pcha_optimal_parameters}} which  conducts a grid search with different values 
 of these parameters and returns the values which minimize the SSE after a fixed number of iterations (\code{testing_iters}, default=10).
}
\section{Evaluate the quality of Archetypal Analysis}{
By using function \code{\link{check_Bmatrix}} we can evaluate the overall quality of 
applied method and algorithm. Quality can be considered high:
\enumerate{
  \item if every archetype is being created by a small number of data points
  \item if relevant weights are not numerically insignificant
}
Of course we must take into account the SSE and VarianceExplained, but if we have to compare two solutions
with similar termination status, then we must choose that of the simplest B matrix form.
} 
\section{Resampling}{
The package includes a function for resampling (\code{\link{grouped_resample}}) which may be used for standard bootstrapping or for subsampling. 
This function allows samples to be drawn with or without replacement, by groups and with or without Dirichlet weights. 
This provides a variety of options for researchers who wish to correct sample biases, estimate empirical confidence intervals, 
and/or subsample large data sets. 
}
\section{Post-run tools}{
Except from \code{\link{check_Bmatrix}} there exist next functions for checking the convergence process itself and
for examining the local neighborhoud of archetypes:
\enumerate{
  \item The function  \code{\link{study_AAconvergence}} analyzes the history of iterations done and produces
  a multi-panel plot showing the steps and quality of the convergence to the final archetypes.
  \item By setting the desired number \code{npoints} as argument in function \code{\link{find_closer_points}} 
  we can then find the data points that are in the local neighborhood of each archetype. 
  This allows us to study the properties of the solution or manually choose an initial approximation to search for a better fit.
}
}
\note{
Bug reports and feature requests can be sent to
  \email{dchristop@econ.uoa.gr} or \cr \email{dem.christop@gmail.com}.
}
\author{
\strong{Maintainer}: Demetris Christopoulos \email{dchristop@econ.uoa.gr}

Other contributors:
\itemize{
  \item David Midgley \email{david.midgley@insead.edu} [contributor]
  \item Sunil Venaik \email{s.venaik@business.uq.edu.au} [contributor]
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