\name{dirichlet_sample}
\alias{dirichlet_sample}
\title{
Function which performs Dirichlet sampling
}
\description{
It uses Dirichlet weights for creating sub-samples of initial data set.
}
\usage{
dirichlet_sample(in_data = NULL, sample_size = NULL,
replacement = NULL, rseed = NULL)
}
\arguments{
  \item{in_data}{
The initial data frame that must be re-sampled. It must contain:
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
  \item{rseed}{
The random seed that will be used for setting initial A matrix. Useful for reproducible results
}
}
\value{
It returns a data frame with exactly the same variables as the initial one, 
except that group variable has now only the given value from input data frame.
}
\author{
David Midgley
}
\seealso{
\code{\link{grouped_resample}}
}
\examples{
## Load absolute temperature data set:
data("AbsoluteTemperature")
df=AbsoluteTemperature
## Find portions for climate zones
pcs=table(df$z)/dim(df)[1]
## Choose the approximate size of the new sample and compute resample sizes
N=1000
resamplesizes=as.integer(round(N*pcs))
sum(resamplesizes)
## Create the grouping matrix
groupmat=data.frame("Group_ID"=1:4,"Resample_Size"=resamplesizes)
groupmat
## Dirichlet resampling:
resample_dirichlet <- grouped_resample(in_data = df,grp_vector = "z",
                      grp_matrix = groupmat,replace = FALSE,
                      option = "Dirichlet", rseed = 20191220)
cat(dim(resample_dirichlet),"\n")
}
\keyword{resampling}
\keyword{Dirichlet}
