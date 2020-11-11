\name{study_AAconvergence}
\alias{study_AAconvergence}
\title{Function which studies the convergence of
Archetypal Analysis when using the PCHA algorithm}
\description{
First it finds an AA solution under given arguments while storing 
all iteration history (\code{save_history = TRUE}).
Then it computes the LOWESS [1] of SSE and its relevant UIK point [2].
Study is performed for iterations after that point.
The list of B-matrices and archetypes that were found are stored.
The archetypes are being aligned, while the B-matrices
are used for computing the used rows-weights,
leading rows-weights and maybe percentage of used rows on Convex Hull.
The Aitken SSE extrapolation plus the relevant error are computed.
The order and rate of convergence  are estimated.
Finally a multi-plot panel is being created if asked.
}
\usage{
study_AAconvergence(df, kappas, method = "projected_convexhull", 
                    rseed = NULL, chvertices = NULL, plot = TRUE, ...)
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

\item{rseed}{The random seed that will be used for setting initial A matrix. Useful for reproducible results.}

\item{chvertices}{The vector of rows which represents the vertices for Convex Hull (if available)}

\item{plot}{If it is TRUE, then a panel of useful plots is created}

\item{...}{Other arguments to be passed to function \code{\link{archetypal}}, except \code{save_history} 
which must always be TRUE}
}
\details{
If we take natural logarithms at the next approximate equation
\deqn{\epsilon_{n+1} = c\epsilon_{n}^p}{%
        e_{n+1} = c e_{n}^p}
for \eqn{n = 1, 2, 3, \ldots}, then we'll find
\deqn{\log(\epsilon_{n+1}) = \log(c)+p\log(\epsilon_{n})}{%
        log(e_{n+1}) = log(c) + p log(e_{n})}
Thus a reasonable strategy for estimating order p and rate c is to perform a linear regression
on above errors, after a selected iteration. 
That is the output of \code{order_estimation} and \code{rate_estimation}.
}
\value{
A list with members:
\enumerate{
\item SSE, a vector of all SSE from all AA iterations
\item SSE_lowess, a vector of LOWESS values for SSE 
\item UIK_lowess, the UIK point [2] of SSE_lowess 
\item aitken, a data frame of Aitken [3] extrapolation and error for SSE after UIK_lowess iteration 
\item order_estimation, the last term in estimating order of convergence, page 56 of [4], 
by using SSE after UIK_lowess iteration 
\item rate_estimation, the last term in estimating rate of convergence, page 56 of [4], 
by using SSE after UIK_lowess iteration 
\item significance_estimations, a data frame with standard errors and statistical significance for estimations
\item used_on_convexhull, the \% of used rows which lie on Convex Hull (if given),
as a sequence for iterations after UIK_lowess one
\item aligned_archetypes, the archetypes after UIK_lowess iteration are being aligned
by using \code{\link{align_archetypes_from_list}}. The history of archetypes creation.
\item solution_used, the AA output that has been used. Some times useful, especially for big data.
}
}
\examples{
{\donttest{
# Load data "wd2"
data(wd2)
ch = chull(wd2)
sa = study_AAconvergence(df = wd2, kappas = 3, rseed = 20191119,
                         verbose = FALSE, chvertices = ch)
names(sa)
# [1] "SSE"                      "SSE_lowess"               "UIK_lowess"              
# [4] "aitken"                   "order_estimation"         "rate_estimation"         
# [7] "significance_estimations" "used_on_convexhull"       "aligned_archetypes"      
# [10] "solution_used"        
# sse=sa$SSE
# ssel=sa$SSE_lowess
sa$UIK_lowess
# [1] 36
# sa$aitken
sa$order_estimation
# [1] 1.007674
sa$rate_estimation
# [1] 0.8277613
sa$significance_estimations
#        estimation   std.error   t.value      p.value
# log(c) -0.1890305 0.014658947 -12.89523 5.189172e-12
# p       1.0076743 0.001616482 623.37475 3.951042e-50
# sa$used_on_convexhull
# sa$aligned_archetypes
data.frame(sa$solution_used[c("SSE","varexpl","iterations","time")])
#        SSE   varexpl iterations time
# 1 1.717538 0.9993186         62 8.39
# Close plot device:
dev.off()
}
}
}
\references{
[1] Cleveland, W. S. (1979) Robust locally weighted regression and smoothing scatterplots. J. Amer. Statist. Assoc. 74, 829--836.\cr

[2] Christopoulos, Demetris T., Introducing Unit Invariant Knee (UIK) As an Objective Choice for
 Elbow Point in Multivariate Data Analysis Techniques (March 1, 2016). 
Available at SSRN: http://dx.doi.org/10.2139/ssrn.3043076 \cr

[3] Aitken, A. "On Bernoulli's numerical solution of algebraic equations", Proceedings of the Royal Society of Edinburgh (1926) 46 pp. 289-305. \cr

[4] Atkinson, K. E.,An Introduction to Numerical Analysis, Wiley & Sons,1989 
}
\seealso{
\code{\link{check_Bmatrix}}
}