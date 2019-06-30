\name{dirichlet_sample}
\alias{dirichlet_sample}
\title{
Function which performs Dirichlet sampling
}
\description{
It uses Dirichlet weights for creating sub-smaples of initial data set.
}
\usage{
dirichlet_sample(in_data = NULL, sample_size = NULL, replacement = NULL)
}
\arguments{
  \item{in_data}{
The initial data frame that must be re-smapled. It must contain:
 \enumerate{
  \item an ID variable
  \item the variables of interest
  \item a grouping variable
 }
}
  \item{sample_size}{
An integer for the size of the new sample
}
  \item{replacement}{
A logical input: TRUE/FALSE if replacement should be used or not, respectively
}
}
\details{
...
}
\value{
It returns a data frame with exactly the same variables as the initial one, 
except that group variable has now only the given value from input data frame.
}
\references{
...
}
\author{
David F. Midgley
}
\note{
...
}
\seealso{
\code{\link{grouped_resample}}
}
\examples{
data("example_data")
data("test_matrix")
resample_dirichlet <- grouped_resample(in_data = example_data,grp_vector = "group",
                      grp_matrix = test_matrix,replace = FALSE,option = "Dirichlet")
cat(dim(resample_dirichlet),"\n")
}
\keyword{resampling}
\keyword{Dirichlet}
