\name{archetypal}
\alias{archetypal}
\title{archetypal: Finds the archetypal analysis of a data frame by using a modified version of PCHA algorithm}
\description{
Performs archetypal analysis by using Convex Hull approximation under a full control of all algorithmic parameters. 
}
\usage{
archetypal(df, kappas, initialrows = NULL,
  method = "projected_convexhull", nprojected = 2, npartition = 10,
  nfurthest = 10, maxiter = 2000, conv_crit = 1e-06,
  var_crit = 0.9999, verbose = TRUE, rseed = NULL, aupdate1 = 25,
  aupdate2 = 10, bupdate = 10, muAup = 1.2, muAdown = 0.5,
  muBup = 1.2, muBdown = 0.5, SSE_A_conv = 1e-09,
  SSE_B_conv = 1e-09, save_history = FALSE, nworkers = NULL)
}
\arguments{
\item{df}{The data frame with dimensions n x d}

\item{kappas}{The number of archetypes}

\item{initialrows}{The initial set of rows from data frame that will be used for starting algorithm}

\item{method}{The method that will be used for computing initial approximation:
\enumerate{
  \item projected_convexhull, see  \code{\link{find_outmost_projected_convexhull_points}} 
  \item convexhull, see \code{\link{find_outmost_convexhull_points}} 
  \item partitioned_convexhull, see \code{\link{find_outmost_partitioned_convexhull_points}} 
  \item furthestsum, see \code{\link{find_furthestsum_points}} 
  \item outmost, see \code{\link{find_outmost_points}} 
  \item random, a random set of kappas points will be used
}}

\item{nprojected}{The dimension of the projected subspace for \code{\link{find_outmost_projected_convexhull_points}}}

\item{npartition}{The number of partitions for \code{\link{find_outmost_partitioned_convexhull_points}}}

\item{nfurthest}{The number of times that \code{\link{FurthestSum}}  algorithm will be applied}

\item{maxiter}{The maximum number of iterations for main algorithm application}

\item{conv_crit}{The SSE convergence criterion of termination: iterate until |dSSE|/SSE<conv_crit}

\item{var_crit}{The Variance Explained (VarExpl) convergence criterion of termination: iterate until VarExpl<var_crit}

\item{verbose}{If it is set to TRUE, then both initialization and iteration details are printed out}

\item{rseed}{The random seed that will be used for seting initial A matrix. Useful for reproducible results.}

\item{aupdate1}{The number of intial applications of Aupdate for improving the initially randomly selected A matrix}

\item{aupdate2}{The number of Aupdate applications in main iteration}

\item{bupdate}{The number of Bupdate applications in main iteration}

\item{muAup}{The factor (>1) by which muA is multiplied when it holds SSE<=SSE_old(1+SSE_A_conv)}

\item{muAdown}{The factor (<1) by which muA is multiplied when it holds SSE>SSE_old(1+SSE_A_conv)}

\item{muBup}{The factor (>1) by which muB is multiplied when it holds SSE<=SSE_old(1+SSE_B_conv)}

\item{muBdown}{The factor (<1) by which muB is multiplied when it holds SSE>SSE_old(1+SSE_B_conv)}

\item{SSE_A_conv}{The convergence value used in SSE<=SSE_old(1+SSE_A_conv). 
Warning: there exists a Matlab crash sometimes after setting this to 1E-16 or lower}

\item{SSE_B_conv}{The convergence value used in SSE<=SSE_old(1+SSE_A_conv).
Warning: there exists a Matlab crash sometimes after setting this to 1E-16 or lower}

\item{save_history}{If set TRUE, then iteration history is being saved for further use}

\item{nworkers}{The number of logical processors that will be used for
parallel computing (usually it is the double of available physical cores)}
}
\value{
A list with members:
\enumerate{
\item BY, the \eqn{kappas \times d} matrix of archetypes found
\item A, the \eqn{n \times kappas} matrix such that Y ~ Y - ABY or Frobenius norm ||Y-ABY|| is minimum
\item B, the \eqn{kappas \times n} matrix such that Y ~ Y - ABY or Frobenius norm ||Y-ABY|| is minimum
\item SSE, the sum of squared error SSE = ||Y-ABY||^2
\item varexpl, the Variance Explained =(SST-SSE)/SST where SST is the total sum of squares for data set matrix
\item initialsolution, the initially used set of rows from data frame in order to start the algorithm
\item freqstable, the frequency table for all found rows, if it is available.
\item iterations, the number of main iterations done by algorithm
\item time, the time in seconds that was spent from entire run 
\item converges, if it is TRUE, then convegence was achieved before the end of maximum allowed iterations
\item nAup, the total number of times when it was SSE<=SSE_old(1+SSE_A_conv) in Aupdate processes. Useful for debugging purposes.
\item nAdown The total number of times when it was SSE>SSE_old(1+SSE_A_conv) in Aupdate processes. Useful for debugging purposes.
\item nBup, the total number of times when it was SSE<=SSE_old(1+SSE_B_conv) in Bupdate processes. Useful for debugging purposes.
\item nBdown, the total number of times when it was SSE>SSE_old(1+SSE_A_conv in Bupdate processes. Useful for debugging purposes.
\item run_results, a list of iteration related details: SSE, varexpl, time, B, BY for all iterations done.
}
}
\examples{
# Create a small 2D data set from 3 corner-points:
p1=c(1,2);p2=c(3,5);p3=c(7,3) 
dp=rbind(p1,p2,p3);dp
set.seed(916070)
pts=t(sapply(1:20, function(i,dp){
  cc=runif(3)
  cc=cc/sum(cc)
  colSums(dp*cc)
},dp))
df=data.frame(pts)
colnames(df)=c("x","y")
# Run AA:
aa=archetypal(df=df,kappas = 3,verbose=FALSE,save_history = TRUE)
# Archetypes:
archs=data.frame(aa$BY)
archs
# See main results:
names(aa)
aa[c("SSE","varexpl","iterations","time")]
# See history of iterations:
names(aa$run_results)
}
\references{

[1] M Morup and LK Hansen, "Archetypal analysis for machine learning and data mining", Neurocomputing (Elsevier, 2012). https://doi.org/10.1016/j.neucom.2011.06.033.

[2] Source: http://www.mortenmorup.dk/index_files/Page327.htm , last accessed 2019-06-07

}
