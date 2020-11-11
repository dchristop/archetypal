find_furthestsum_points=function(df, kappas, nfurthest = 100, nworkers = NULL,
                                 sortrows = TRUE, doparallel = TRUE){
  # Check logical processors and define "nworkers" if necessary
  if(doparallel){
    nwa=parallel::detectCores()
    if(is.null(nworkers)){
      ifelse(nwa>2,{nworkers=nwa-2},{nworkers=nwa})
    }else{
      if(nworkers>nwa){nworkers=nwa}
    }
  }
  # Initialize
  Y=as.matrix(df)
  NN=1:dim(Y)[1]
  # Do 'nfurthest' applications of FurthestSum algorithm 
  if(kappas!=1){
    #Function:
    ifelse(sortrows,
           {runfs=function(j,Y,kappas,NN){sort(FurthestSum(Y=Y,kappas=kappas,irows=sample(NN,1)))}},
           {runfs=function(j,Y,kappas,NN){FurthestSum(Y=Y,kappas=kappas,irows=sample(NN,1))}})
    environment(runfs) <- .GlobalEnv
    # Choose computation mode and proceed
    if(doparallel){
      #In parallel
      cl <- makeCluster(nworkers);registerDoParallel(cl)
      clusterEvalQ(cl=cl, library("Matrix"))
      clusterExport(cl=cl,varlist='FurthestSum')
      tfsum=t(parallel::parSapply(cl=cl,1:nfurthest,runfs,Y=Y,kappas=kappas,NN=NN))
      stopCluster(cl)
      #
    }else{
      #In serial
      tfsum=t(sapply(1:nfurthest,runfs,Y=Y,kappas=kappas,NN=NN))
      #
    }
    # Find outermost points
    outmostrows=do.call(c,as.list(tfsum))
    di=as.data.frame(table(outmostrows),stringsAsFactors = F)
    di=di[order(-di$Freq),]
    di$outmostrows=as.integer(di$outmostrows)
    di$FreqPerCent=di$Freq/sum(di$Freq)
    di$CumFreqPerCent=cumsum(di$FreqPerCent)
    rownames(di)=1:dim(di)[1]
    #
    rowsmaxall=di$outmostrows
    #
    rowsmaxkappas=rowsmaxall[1:kappas]
    #
    out=list("outmost"=rowsmaxkappas,"outmostall"=rowsmaxall,"outmostfrequency"=di)
  }else{
    cat('For number of archetypes k = 1 the FurthestSum algorithm always gives the first point of data set','\n')
    freqstable=data.frame("outmostrows"=1,"Freq"=1, "FreqPerCent"=1,"CumFreqPerCent"=1)
    out=list("outmost"=1L,"outmostall"=1L,"outmostfrequency"=freqstable)
  }
  #OK all
  #Return list of results
  return(out)
}