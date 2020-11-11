grouped_resample <- function (in_data = NULL, grp_vector = NULL, grp_matrix = NULL, 
                             replace = FALSE, option = "Simple", number_samples = 1,
                             nworkers = NULL, rseed = NULL){
  #
  ######################################################
  # Set seed for reproducible results or set it randomly
  ######################################################
  #
  if(!is.null(rseed)){
    rseed=rseed
    set.seed(rseed)
  }else{
    rseed=sample(1:10^9,1)
  }
  #
  ##################################################################
  # Define number of workers if necessary or fix a wrong given value
  ##################################################################
  #
  nwall=parallel::detectCores()
  if(is.null(nworkers)){
    if(nwall<=2){nworkers=nwall}else{nworkers=nwall-2}
  }else{
    if(nworkers>nwall){nworkers=nwall}
  }
  #
  #########################################
  # Work in parallel if number_smaples > 1
  #########################################
  #
  resamples=list()
  #
  if(number_samples==1){
    for (ii in 1:nrow(grp_matrix)) {
      a_group <- in_data[which(in_data[, grp_vector] == ii), 
                         ]
      grp_sample_size <- grp_matrix[which(grp_matrix[, "Group_ID"] == 
                                            ii), "Resample_Size"]
      if (option == "Simple"){
        if(!is.null(rseed)){set.seed(rseed+ii)}
        a_grp_sample <- a_group[sample(1:nrow(a_group), size = grp_sample_size, replace = replace), ]
      }
      else{
        a_grp_sample <- dirichlet_sample(in_data = a_group, sample_size = grp_sample_size,
                                         replacement = replace, rseed = rseed+ii)
      }
      ifelse(ii == 1, {a_resample <- a_grp_sample}, {a_resample <- rbind(a_resample, a_grp_sample)})
      resamples[[1]]=a_resample
    }
  }else if(number_samples>1){
    # Function:
    runresampling=function(i, in_data, grp_vector, grp_matrix, replace, option, rseed){
      for (ii in 1:nrow(grp_matrix)) {
        a_group <- in_data[which(in_data[, grp_vector] == ii), 
                           ]
        grp_sample_size <- grp_matrix[which(grp_matrix[, "Group_ID"] == ii), "Resample_Size"]
        if (option == "Simple"){
          if(!is.null(rseed)){set.seed(rseed+i)}
          a_grp_sample <- a_group[sample(1:nrow(a_group), size = grp_sample_size, replace = replace), ]
        }
        else {
          a_grp_sample <- dirichlet_sample(in_data = a_group, sample_size = grp_sample_size, 
                                           replacement = replace, rseed = rseed+i)
        }
        if (ii == 1){a_resample <- a_grp_sample}else{a_resample <- rbind(a_resample, a_grp_sample)}
      }
      return(a_resample)
    }
    environment(runresampling) <- .GlobalEnv
    tp1=Sys.time()
    cl <- makeCluster(nworkers);registerDoParallel(cl)
    # clusterSetRNGStream(cl = cl, rseed)
    clusterEvalQ(cl=cl,list(library("archetypal")))
    resamples=parallel::parLapply(cl=cl, 1:number_samples, runresampling,
                                  in_data, grp_vector, grp_matrix, replace, option, rseed)
    stopCluster(cl)
    tp2=Sys.time();cat(" ","\n");print(tp2-tp1) 
    #
  }else{
    stop("The number of samples must be a positive integer")
  }
  # Return the list of re-samples
  return(resamples)
}