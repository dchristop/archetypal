check_Bmatrix=function(B,print.details=TRUE){
  # B from next AA approximation:
  # Y ~ A B Y
  # The archetypes are the matrix B Y
  pw=list()
  ww=list()
  for(i in 1:dim(B)[1]){
    if(print.details){
      cat(paste0("Archetype ",i," is a mixture of only next rows with weights:"),"\n")
      cat(" ","\n")
    }
    dbi=B[i,]
    names(dbi)=1:dim(B)[2]
    rni=dbi[order(dbi,decreasing = T)]
    rni=rni[rni!=0]
    pw[[i]]=as.integer(names(rni))
    ww[[i]]=rni
    if(print.details){
      print(rni)
      cat("Weights add to:","\n")
      print(sum(rni))
      cat(" ","\n")
      cat(" ","\n")
    }
  }
  #
  return(list("used_rows"=pw,"used_weights"=ww))
}