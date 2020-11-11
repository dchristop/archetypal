dirichlet_sample <- function(in_data=NULL,sample_size=NULL,replacement=NULL, rseed = NULL) 
{
  ## function for dirichlet sampling 
  #
  ###################################
  # Set seed for reproducible results
  ###################################
  #
  if(!is.null(rseed)){rseed=rseed}
  #
  # Assign dirichlet weights
  d_weights <- matrix(rexp(nrow(in_data),1),ncol=nrow(in_data),byrow = TRUE)
  d_weights <- d_weights/rowSums(d_weights)
  # care here with "x"
  x <- 1:nrow(in_data)
  # Sample
  if(!is.null(rseed)){set.seed(rseed)}
  a_sample <- in_data[sample(x = x,size = sample_size, replace = replacement, prob = d_weights[1,]),]
  cat("dirichlet sample of",sample_size,"cases drawn\n")
  return(a_sample)
}