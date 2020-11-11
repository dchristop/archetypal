\name{find_pcha_optimal_parameters}
\alias{find_pcha_optimal_parameters}
\title{Finds the optimal updating parameters to be used for the PCHA algorithm}
\description{
After creating a grid on the space of (mu_up, mu_down) it runs \code{\link{archetypal}} by using a given \code{method} & other running options passed by ellipsis (...) and finally finds those values which minimize the SSE at the end of \code{testing_iters} iterations (default=10).
}
\usage{
find_pcha_optimal_parameters(df, kappas, method = "projected_convexhull", 
testing_iters = 10, nworkers = NULL, nprojected = 2, npartition = 10,
nfurthest = 100, sortrows = FALSE,
mup1 = 1.1, mup2 = 2.50, mdown1 = 0.1, mdown2 = 0.5, nmup = 10, nmdown = 10,
rseed = NULL, plot = FALSE, ...)
}
\arguments{
\item{df}{The data frame with dimensions n x d}

\item{kappas}{The number of archetypes}

\item{method}{The method that will be used for computing initial approximation:
\enumerate{
  \item projected_convexhull, see  \code{\link{find_outmost_projected_convexhull_points}} 
  \item convexhull, see \code{\link{find_outmost_convexhull_points}} 
  \item partitioned_convexhull, see \code{\link{find_outmost_partitioned_convexhull_points}} 
  \item furthestsum, see \code{\link{find_furthestsum_points}} 
  \item outmost, see \code{\link{find_outmost_points}} 
  \item random, a random set of kappas points will be used
}}

\item{testing_iters}{The maximum number of iterations to run for every pair (mu_up, mu_down) of parameters}

\item{nworkers}{The number of logical processors that will be used for
parallel computing (usually it is the double of available physical cores)}

\item{nprojected}{The dimension of the projected subspace for \code{\link{find_outmost_projected_convexhull_points}}}

\item{npartition}{The number of partitions for \code{\link{find_outmost_partitioned_convexhull_points}}}

\item{nfurthest}{The number of times that \code{\link{FurthestSum}}  algorithm will be applied}

\item{sortrows}{If it is TRUE, then rows will be sorted in \code{\link{find_furthestsum_points}}}

\item{mup1}{The minimum value of mu_up, default is 1.1}

\item{mup2}{The maximum value of mu_up, default is 2.5}

\item{mdown1}{The minimum value of mu_down, default is 0.1}

\item{mdown2}{The maximum value of mu_down, default is 0.5}

\item{nmup}{The number of points to be taken for [mup1,mup2], default is 10}

\item{nmdown}{The number of points to be taken for [mdown1,mdown2]}

\item{rseed}{The random seed that will be used for setting initial A matrix. Useful for reproducible results}

\item{plot}{If it is TRUE, then a 3D plot for (mu_up, mu_down, SSE) is created}

\item{...}{Other arguments to be passed to function \code{\link{archetypal}}}

}
\value{
A list with members:
\enumerate{
\item mu_up_opt, the optimal found value for muAup and muBup
\item mu_down_opt, the optimal found value for muAdown and muBdown
\item min_sse, the minimum SSE which corresponds to (mu_up_opt,mu_down_opt)
\item seed_used, the used random seed, absolutely necessary for reproducing optimal results
\item method_used, the method that was used for creating the initial solution
\item sol_initial, the initial solution that was used for all grid computations
\item testing_iters, the maximum number of iterations done by every grid computation
}
}
\examples{
{\donttest{
data("wd25")
out = find_pcha_optimal_parameters(df = wd25, kappas = 5, rseed = 2020)
# Time difference of 30.91101 secs
# mu_up_opt mu_down_opt     min_sse 
# 2.188889    0.100000    4.490980  
# Run now given the above optimal found parameters:
aa = archetypal(df = wd25, kappas = 5,
                initialrows = out$sol_initial, rseed = out$seed_used,
                muAup = out$mu_up_opt, muAdown = out$mu_down_opt,
                muBup = out$mu_up_opt, muBdown = out$mu_down_opt)
aa[c("SSE", "varexpl", "iterations", "time" )]
# $SSE
# [1] 3.629542
# 
# $varexpl
# [1] 0.9998924
# 
# $iterations
# [1] 146
# 
# $time
# [1] 21.96
# Compare it with a simple solution (time may vary)
aa2 = archetypal(df = wd25, kappas = 5, rseed = 2020)
aa2[c("SSE", "varexpl", "iterations", "time" )]
# $SSE
# [1] 3.629503
# 
# $varexpl
# [1] 0.9998924
# 
# $iterations
# [1] 164
# 
# $time
# [1] 23.55
## Of course the above was a "toy example", if your data has thousands or million rows,
## then the time reduction is much more conspicuous.
# Close plot device:
dev.off()
}
}
}
\seealso{
\code{\link{find_closer_points}}
}