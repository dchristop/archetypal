\name{find_closer_points}
\alias{find_closer_points}
\title{Function which finds the data points that are closer to the archetypes during all
iterations of the algorithm PCHA}
\description{
This function runs the PCHA algorithm and finds the data points that are in the local neighborhood of each archetype. The size of the neighborhood is user defined (\code{npoints}). This allows us to study the properties of the solution or manually choose an initial approximation to search for a better fit. 
}
\usage{
find_closer_points(df, kappas, usedata = FALSE, npoints = 2, 
                     nworkers = NULL, rseed = NULL, 
                     verbose = FALSE, doparallel = FALSE, ...)
}
\arguments{
\item{df}{The data frame with dimensions n x d}

\item{kappas}{The number of archetypes}

\item{usedata}{If it is TRUE, then entire data frame will be used, if \code{doparallel = TRUE}}

\item{npoints}{The number of closer points to be estimated}

\item{nworkers}{The number of logical processors that will be used, if \code{doparallel = TRUE}}

\item{rseed}{The random seed that will be used for random generator. Useful for reproducible results.}

\item{verbose}{If it is set to TRUE, then details will be printed, except from \code{archetypal}}

\item{doparallel}{If it is set to TRUE, then parallel processing will be performed}

\item{...}{Other arguments to be passed to \code{archetypal} except internally used \code{save_history = TRUE}
and \code{verbose = FALSE}. This is essential for using optimal parameters found by \code{find_pcha_optimal_parameters}}
}
\value{
A list with members:
\enumerate{
  \item rows_history, a list with \code{npoints} rows used that are closer to each archetype
  for each iteration done by algorithm
  \item iter_terminal, iteration after which rows closer to archetypes do not change any more
  \item rows_closer, the rows closer to archetypes by means of Euclidean distance and are fixed
  after \code{iter_terminal} iteration
  \item rows_closer_matrix, a matrix with ]code{npoints} rows  which are closer to each archetype
  \item solution_used, the AA output that has been used. Some times useful, especially for big data.
}
}
\examples{
{\donttest{
# Load data "wd2"
data("wd2")
yy = find_closer_points(df = wd2, kappas = 3, npoints = 2, nworkers = 2)
yy$rows_history
yy$iter_terminal
yy$rows_closer
yy$rows_closer_matrix
yy$solution_used$BY
}
}
}
\seealso{
\code{\link{check_Bmatrix}}, \code{\link{study_AAconvergence}}
}