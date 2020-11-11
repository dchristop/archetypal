find_optimal_kappas=function(df, maxkappas = 15, method = "projected_convexhull",
                             ntrials = 10, nworkers = NULL , ...){
  # Fix number of maxkappas if necessary
  if(maxkappas<6){maxkappas=6}
  # Define number of workers if necessary
  #
  if(is.null(nworkers)){
    nwall=parallel::detectCores()
    if(nwall<=2){nworkers=nwall}else{nworkers=nwall-2}
  }
  #
  #################################################################
  # Run for kappas=1:maxkappas and 'ntrials' times for each kappas:
  #################################################################
  # Function ...
  runaaN=function(j, k, df, method, ...){
    archetypal(df = df, kappas = k, method = method, ...)
  }
  environment(runaaN) <- .GlobalEnv
  # Cases ...
  if(ntrials>1){
    # Run 'ntrials' times for each kappas, in parallel
    #
    tp1=Sys.time()
    dkappas=list()
    for(kappas in 1:maxkappas){
      cat(paste0(kappas,' '))
      cl <- makeCluster(nworkers);registerDoParallel(cl);
      clusterEvalQ(cl=cl,list(library("archetypal")))
      dkappas[[kappas]]=parallel::parLapply(cl = cl, 1:ntrials, runaaN,
                        k = kappas, df = df, method = method, ...)
      stopCluster(cl)
    }
    tp2=Sys.time();cat(" ","\n");print(tp2-tp1) 
    names(dkappas)=paste0("kappas_",1:maxkappas)
    #
  }else{
    # Run only once for each kappas but in parallel for all of them
    #
    tp1=Sys.time()
    dkappas=list()
    cl <- makeCluster(nworkers);registerDoParallel(cl);
    clusterEvalQ(cl=cl,list(library("archetypal")))
    dkappas=parallel::parLapply(cl=cl,1:maxkappas,function(k,df,method, ...){
      archetypal(df = df, kappas = k, method = method, ...)},
      df = df, method = method, ...)
    stopCluster(cl)
    tp2=Sys.time();cat(" ","\n");print(tp2-tp1) 
    names(dkappas)=paste0("kappas_",1:maxkappas)
    # 
  }
  # SSE and VarExpl
  if(ntrials>1){
    # SSE
    msse=data.frame(t(sapply(dkappas, function(x){sapply(x, function(y){y$SSE})})))
    colnames(msse)=1:dim(msse)[2]
    msse$kappas=1:(dim(msse)[1])
    # VarExpl
    mvex=data.frame(t(sapply(dkappas, function(x){sapply(x, function(y){y$varexpl})})))
    colnames(mvex)=1:dim(mvex)[2]
    mvex$kappas=1:(dim(mvex)[1])
  }else{
    msse=t(sapply(dkappas, function(x){x$SSE}))
    mvex=t(sapply(dkappas, function(x){x$varexpl}))
  }
  #
  ###################################
  # Find UIK FOR ALL AVAILABLE VALUES
  ###################################
  # Cases ...
  if(ntrials>1){
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
  }else{
    xsse=1:maxkappas;ysse=msse[1,];
    dsse=data.frame("k"=xsse,"sse"=ysse)
    dsse1=dsse;dsse1$sse=100*dsse$sse/dsse$sse[1];
    dbest=dsse;dbest1=dsse1
    ysse1=dsse1$sse
    xbest=xsse;ybest=ysse;
    ybest1=ysse1
  }
  #
  ################################
  # Find UIK for all possible SSEs
  ################################
  #
  knee1=uik(xsse,ysse)
  knee2=uik(xsse,ysse1)
  knee3=uik(xbest,ybest)
  knee4=uik(xbest,ybest1)
  #
  ##############################
  # FIND D2UIK FOR BEST FIT ONLY
  ##############################
  #
  uik2=d2uik(xbest,ybest)
  #
  ###################################################
  # Return all. Set optimal from best fit solution
  ###################################################
  #
  out=list("all_sse"=dsse,"all_sse1"=dsse1,"bestfit_sse"=dbest,"bestfit_sse1"=dbest1,
           "all_kappas"=c(knee1,knee2,knee3,knee4),"d2uik"=uik2,"optimal_kappas"= knee3)
  return(out)
  #
}
