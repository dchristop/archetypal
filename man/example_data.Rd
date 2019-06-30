\name{example_data}
\alias{example_data}
\docType{data}
\title{
An example data set for resampling functions
}
\description{
A 100K x 12 data frame for demonstrating resampling functions \code{\link{dirichlet_sample}} and \code{\link{grouped_resample}}.
}
\usage{data("example_data")}
\format{
  A data frame with 100000 observations on the following 12 variables.
  \describe{
    \item{\code{id}}{a numeric vector with the ID variable}
    \item{\code{opn01}}{a numeric vector}
    \item{\code{opn02}}{a numeric vector}
    \item{\code{opn03}}{a numeric vector}
    \item{\code{opn04}}{a numeric vector}
    \item{\code{opn05}}{a numeric vector}
    \item{\code{opn06}}{a numeric vector}
    \item{\code{opn07r}}{a numeric vector}
    \item{\code{opn08}}{a numeric vector}
    \item{\code{opn09r}}{a numeric vector}
    \item{\code{opn10}}{a numeric vector}
    \item{\code{group}}{a numeric vector with the grouping variable}
  }
}
\details{
...
}
\source{
...
}
\references{
...
}
\examples{
data(example_data)
table(example_data$group)
}
\keyword{datasets}
