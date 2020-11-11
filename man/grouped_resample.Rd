\name{grouped_resample}
\alias{grouped_resample}
\title{
Function for performing simple or Dirichlet resampling
}
\description{
The function may be used for standard bootstrapping or for subsampling, see [1]. 
This function allows samples to be drawn with or without replacement, by groups and with or without 
Dirichlet weights, see [2]. This provides a variety of options for researchers who wish
to correct sample biases, estimate empirical confidence intervals, and/or subsample large data sets. 
}
\usage{
grouped_resample(in_data = NULL, grp_vector = NULL, grp_matrix = NULL, 
                 replace = FALSE, option = "Simple", number_samples = 1, 
                 nworkers = NULL, rseed = NULL)
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
  \item{grp_vector}{
The grouping variable of the data frame, defined under the name 'group' for example
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
  \item{number_samples}{
The number of samples to be created. If it is greater than one, then parallel processing is used.
}
  \item{nworkers}{
The number of logical processors that will be used for parallel computing
(usually it is the double of available physical cores)
}
  \item{rseed}{
The random seed that will be used for sampling. Useful for reproducible results
}
}
\value{
It returns a list of \code{mumber_samples} data frames with exactly the same 
variables as the initial one, except that group variable has now only the given 
value from input data frame.
}
\references{
[1] D. N. Politis, J. P. Romano, M. Wolf, Subsampling (Springer-Verlag, New York, 1999).

[2] Baath R (2018). bayesboot: An Implementation of Rubin's (1981) Bayesian Bootstrap. R
package version 0.2.2, URL https://CRAN.R-project.org/package=bayesboot
}
\author{
David Midgley
}
\seealso{
\code{\link{dirichlet_sample}}
}
\examples{
## Load absolute temperature data set:
data("AbsoluteTemperature")
df <- AbsoluteTemperature
## Find portions for climate zones
pcs <- table(df$z)/dim(df)[1]
## Choose the approximate size of the new sample and compute resample sizes
N <- round(sqrt(nrow(AbsoluteTemperature)))
resamplesizes=as.integer(round(N*pcs))
sum(resamplesizes)
## Create the grouping matrix
groupmat <- data.frame("Group_ID"=1:4,"Resample_Size"=resamplesizes)
groupmat
## Simple resampling:
resample_simple <- grouped_resample(in_data = df, grp_vector = "z",
                                    grp_matrix = groupmat, replace = FALSE, option = "Simple",
                                    number_samples = 1, nworkers = NULL, rseed = 20191220)
cat(dim(resample_simple[[1]]),"\n")
## Dirichlet resampling:
resample_dirichlet <- grouped_resample(in_data = df, grp_vector = "z",
                                       grp_matrix = groupmat, replace = FALSE, option = "Dirichlet",
                                       number_samples = 1, nworkers = NULL, rseed = 20191220)
cat(dim(resample_dirichlet[[1]]),"\n")
##
# ## Work in parallel and create many samples
# ## Choose a random seed
# nseed <- 20191119
# ## Simple
# reslist1 <- grouped_resample(in_data = df, grp_vector = "z", grp_matrix = groupmat,
#                         replace = FALSE, option = "Simple",
#                         number_samples = 10, nworkers = NULL,
#                         rseed = nseed)
# sapply(reslist1, dim)
# ## Dirichlet
# reslist2 <- grouped_resample(in_data = df, grp_vector = "z", grp_matrix = groupmat,
#                         replace = FALSE, option = "Dirichlet",
#                         number_samples = 10, nworkers = NULL,
#                         rseed = nseed)
# sapply(reslist2, dim)
# ## Check for same rows between 1st sample of 'Simple' and 1st sample of 'Dirichlet' ...
# mapply(function(x,y){sum(rownames(x)%in%rownames(y))},reslist1,reslist2)
#
}
\keyword{resampling}
\keyword{Dirichlet}