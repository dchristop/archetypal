\name{test_matrix}
\alias{test_matrix}
\docType{data}
\title{
An example of a grouping matrix
}
\description{
A grouping matrix for demonstrating resampling functions \code{\link{dirichlet_sample}} and \code{\link{grouped_resample}}.
}
\usage{data("test_matrix")}
\format{
  The format is:
 num [1:3, 1:2] 1 2 3 317 698 431
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:2] "Group_ID" "Resample_Size"
}
\details{
...
}
\references{
...
}
\examples{
data(test_matrix)
print(test_matrix)
}
\keyword{datasets}