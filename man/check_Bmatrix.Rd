\name{check_Bmatrix}
\alias{check_Bmatrix}
\title{Function which checks B matrix of Archetypal Analysis  Y ~ A B Y
in order to find the used rows for creating each archetype
and the relevant used weights.}
\description{
Function which checks B matrix of Archetypal Analysis  Y ~ A B Y
in order to find the used rows for creating each archetype
and the relevant used weights.
}
\usage{
check_Bmatrix(B, chvertices = NULL, verbose = TRUE)
}
\arguments{
\item{B}{The \eqn{kappas \times n} matrix such that Y ~ ABY or Frobenius norm ||Y-ABY|| is minimum}
\item{chvertices}{The vector of rows which represent the Convex Hull of data frame}
\item{verbose}{If set to TRUE, then results are printed out.}
}
\value{
A list with members:
\enumerate{
  \item used_rows, a list with used rows for creating each archetype
  \item used_weights, a list with the relevant weights that have been used
  \item leading_rows, the rows for each archetype with greatest weight
  \item leading_weights, the weights of leading rows
  \item used_on_convexhull, the portion of used rows which lie on Convex Hull (if given)
}
}
\examples{
{\donttest{
# Load data "wd2"
data("wd2")
df = wd2
# Run AA:
aa = archetypal(df = df, kappas = 3, verbose = FALSE)
# Check B matrix:
B = aa$B
yy = check_Bmatrix(B, verbose = TRUE)
yy$used_rows
yy$used_weights
yy$leading_rows
yy$leading_weights
# Check if used rows lie on ConvexHull
ch = chull(df)
yy = check_Bmatrix(B, chvertices = ch, verbose = FALSE)
yy$used_on_convexhull
#
}
}
}
\seealso{
\code{\link{archetypal}}, \code{\link{check_Bmatrix}}, \code{\link{find_closer_points}}  \cr
& \code{\link{study_AAconvergence}}
}
