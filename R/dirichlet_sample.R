dirichlet_sample <- function(in_data=NULL,sample_size=NULL,replacement=NULL) 
{
  # function for dirichlet sampling 
  # Assign dirichlet weights
  d_weights <- matrix(rexp(nrow(in_data),1),ncol=nrow(in_data),byrow = TRUE)
  d_weights <- d_weights/rowSums(d_weights)
  # care here with "x"
  x <- 1:nrow(in_data)
  # Sample
  a_sample <- in_data[sample(x = x,size = sample_size, replace = replacement, prob = d_weights[1,]),]
  cat("dirichlet sample of",sample_size,"cases drawn\n")
  return(a_sample=a_sample)
}