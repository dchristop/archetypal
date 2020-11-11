\name{find_optimal_kappas}
\alias{find_optimal_kappas}
\title{Function for finding the optimal number of archetypes}
\usage{
find_optimal_kappas(df, maxkappas = 15, method = "projected_convexhull", 
                    ntrials = 10, nworkers = NULL, ...)
}
\arguments{
\item{df}{The data frame with dimensions \eqn{n \times d}}

\item{maxkappas}{The maximum number of archetypes for which algorithm will be applied}

\item{method}{The method that will be used for computing the initial solution}

\item{ntrials}{The number of times that algorithm will be applied for each kappas}

\item{nworkers}{The number of logical processors that will be used for parallel computing 
 (usually it is the double of available physical cores)}
 
\item{...}{Other arguments to be passed to function \code{\link{archetypal}}}
 
}
\value{
A list with members:
\enumerate{
  \item all_sse, all available SSE for all kappas and all trials per kappas
  \item all_sse1, all available SSE(k)/SSE(1) for all kappas and all trials per kappas
  \item bestfit_sse, only the best fit SSE trial for each kappas
  \item bestfit_sse1, only the best fit SSE(k)/SSE(1) trial for each kappas
  \item all_kappas, the knee point of scree plot for all 4 SSE results
  \item d2uik, the UIK for the absolute values of the estimated best fit SSE second derivatives,
  after using second order forward divided differences approximation
  \item optimal_kappas, the knee point from best fit SSE results
}
}
\description{
Function for finding the optimal number of archetypes in order to apply Archetypal Analysis for a data frame.
}
\details{
After having found the SSE for each kappas, UIK method (see [1]) is used
for estimating the knee or elbow point as the optimal kappas. 
}
\examples{
{\donttest{
# Run may take a while depending on your machine ...
# Load data frame "wd2"
data("wd2")
df = wd2
# Run:
t1 = Sys.time()
yy = find_optimal_kappas(df, maxkappas = 10)
t2 = Sys.time();print(t2-t1)
# Results:
names(yy)
# Best fit SSE:
yy$bestfit_sse 
# Optimal kappas from UIK method:
yy$optimal_kappas
#
}
}
}
\references{
[1] Christopoulos, Demetris T., Introducing Unit Invariant Knee (UIK) As an Objective Choice for
 Elbow Point in Multivariate Data Analysis Techniques (March 1, 2016). 
Available at SSRN: http://dx.doi.org/10.2139/ssrn.3043076
}
\seealso{
\code{\link{archetypal}}
}
