grouped_resample <- function(in_data=NULL,grp_vector=NULL,grp_matrix=NULL,replace=FALSE,
                             option="Simple")
{
  for(ii in 1:nrow(grp_matrix)) {
    a_group <- in_data[which(in_data[,grp_vector]==ii),]
    #cat(ii,dim(a_group),"\n")
    grp_sample_size <- grp_matrix[which(grp_matrix[,"Group_ID"]==ii),"Resample_Size"]
    if(option=="Simple") {
      a_grp_sample <- a_group[sample(1:nrow(a_group),size = grp_sample_size,replace = replace),]
      #cat(dim(a_grp_sample),"\n")
    } else {
      a_grp_sample <- dirichlet_sample(in_data = a_group, sample_size = grp_sample_size, 
                                       replacement = replace)
    }
    if(ii == 1) a_resample <- a_grp_sample else a_resample <- rbind(a_resample,a_grp_sample)
  }
  return(a_resample=a_resample)
}