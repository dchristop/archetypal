find_optimal_kappas=function(df,maxkappas=15,method="projected_convexhull",ntrials=10,nworkers=NULL){
  # Define number of workers if necessary
  if(is.null(nworkers)){
    nwall=parallel::detectCores()
    if(nwall<=2){nworkers=2}else{nworkers=nwall-2}
  }
  # Run for kappas=1:maxkappas and 'ntrials' times for each kappas:
  tp1=Sys.time()
  dkappas=list()
  # Function:
  runaaN=function(j,k,df,method,muoptUP,muoptDOWN){
    archetypal(df=df,kappas=k,method = method, aupdate1=30,aupdate2=5,bupdate=5,
               SSE_B_conv = 0,SSE_A_conv = 0,save_history = FALSE,maxiter = 2000)
  }
  environment(runaaN) <- .GlobalEnv
  # Run 'ntrials' times for each kappas
  for(kappas in 1:maxkappas){
    cat(kappas,'\n')
    cl <- makeCluster(nworkers);registerDoParallel(cl);
    clusterEvalQ(cl=cl,list(library("Matrix"),library("geometry"),library("archetypal")))
    dkappas[[kappas]]=parallel::parLapply(cl=cl,1:ntrials,runaaN,k=kappas,df=df,method=method)
    stopCluster(cl)
  }
  tp2=Sys.time();print(tp2-tp1) 
  names(dkappas)=paste0("kappas_",1:maxkappas)
  # SSE
  msse=data.frame(t(sapply(dkappas, function(x){sapply(x, function(y){y$SSE})})))
  colnames(msse)=1:dim(msse)[2]
  msse$kappas=1:(dim(msse)[1])
  # VarExpl
  mvex=data.frame(t(sapply(dkappas, function(x){sapply(x, function(y){y$varexpl})})))
  colnames(mvex)=1:dim(mvex)[2]
  mvex$kappas=1:(dim(mvex)[1])
  #
  ###########
  # Find UIK 
  ###########
  # ALL AVAILABLE VALUES
  xsse=rep(msse[1,"kappas"],ntrials)
  ysse=as.vector(unlist(msse[1,1:ntrials]))
  for(j in 2:maxkappas){
    xsse=c(xsse,rep(msse[j,"kappas"],ntrials))
    ysse=c(ysse,as.vector(unlist(msse[j,1:ntrials])))
  }
  # Define data frames and variables for SSE
  dsse=data.frame("k"=xsse,"sse"=ysse)
  dsse1=dsse;dsse1$sse=100*dsse$sse/dsse$sse[1];
  ysse1=dsse1$sse
  # Get best fit values
  dbest=aggregate(sse~k,dsse,function(x){imin=which.min(x);c(imin,x[imin])});
  dbest=data.frame("k"=dbest$k,"imin"=as.integer(dbest$sse[,1]),"sse"=dbest$sse[,2])
  xbest=dbest$k
  ybest=dbest$sse
  ybest1=100*ybest/ybest[1]
  dbest1=dbest
  dbest1$sse=ybest1
  # Find UIK for all cases
  knee1=xsse[ede(xsse,ysse,0)[1]]
  knee2=xsse[ede(xsse,ysse1,0)[1]]
  knee3=xbest[ede(xbest,ybest,0)[1]]
  knee4=xbest[ede(xbest,ybest1,0)[1]]
  #
  return(list("all_sse"=dsse,"all_sse1"=dsse1,"bestfit_sse"=dbest,"bestfit_sse1"=dbest1,
              "all_kappas"=c(knee1,knee2,knee3,knee4),"optimal_kappas"= knee3))
  #
}
