\name{grouped_resample}
\alias{grouped_resample}
\title{
Function for performing simple or Dirichlet resampling
}
\description{
This function can create either a simple or a Dirichlet weighted sub-sample, given an initial data frame
and a grouping variable
}
\usage{
grouped_resample(in_data = NULL, grp_vector = NULL, grp_matrix = NULL, 
                 replace = FALSE, option = "Simple")
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
  \item{grp_vector}{
The grouping variable of the data frame, defined under the name 'froup' for example
}
  \item{grp_matrix}{
A matrix that contains
 \enumerate{
  \item the variable 'Group_ID' with entries all the available values of grouping variable 
  \item the variable 'Resample_Size' with the sizes for each sample that will be created per grouping value
 }
}
  \item{replace}{
A logical input: TRUE/FALSE if replacement should be used or not, respectively
}
  \item{option}{
A character input with next possible values
\enumerate{
  \item "Simple", if we want to perform a simple re-sampling
  \item "Dirichlet", if we want to perform a Dirichlet weighted re-sampling
 }
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
\code{\link{dirichlet_sample}}
}
\examples{
## Load example data set and grouping matrix:
data("example_data")
data("test_matrix")
## Simple resampling:
resample_simple <- grouped_resample(in_data = example_data,grp_vector = "group",
                     grp_matrix = test_matrix,replace = FALSE,option = "Simple")
cat(dim(resample_simple),"\n")
## Dirichlet resampling:
resample_dirichlet <- grouped_resample(in_data = example_data,grp_vector = "group",
                      grp_matrix = test_matrix,replace = FALSE,option = "Dirichlet")
cat(dim(resample_dirichlet),"\n")
}
\keyword{resampling}
\keyword{Dirichlet}