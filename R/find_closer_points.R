find_closer_points=function(df, kappas, usedata = FALSE, npoints = 2, 
                              nworkers = NULL, rseed = NULL, 
                              verbose = FALSE, doparallel = FALSE, ...){
  ####################################################
  #Find closer points to archetypes for PCHA algorithm
  ####################################################
  # Convert row names to integers 
  rownames(df)=1:dim(df)[1]
  # Set seed for reproducible results
  if(!is.null(rseed)){rseed=rseed}
  # Detect and set workers if needed
  if(doparallel){
    nwa=detectCores()
    if(is.null(nworkers)){
      ifelse(nwa>2,{nworkers=nwa-2},{nworkers=nwa})
    }else{
      if(nworkers>nwa){nworkers=nwa}
    }
  }
  # Check npoints
  if(npoints>dim(df)[1]){npoints=dim(df)[1]-1}
  # Run AA
  # If variables <=6 then run on CH, otherwise run on initial data frame
  if(dim(df)[2]<=6 * !usedata){
    if(verbose){cat('Work on Convex Hull vertices only, since they were easily computable','\n')}
    ch=unique(do.call(c,as.list(convhulln(df,'FA')$hull)))
    dch=df[ch,]
    aa = archetypal(df = dch, kappas = kappas, rseed = rseed,
                    save_history = TRUE, verbose = FALSE, ...)
    archs=aa$BY
    archslist=aa$run_results$archslist
  }else{
    if(verbose){cat('Work on entire data frame , since Convex Hull was not easily computable','\n')}
    aa = archetypal(df = df, kappas = kappas, rseed = rseed, 
                    save_history = TRUE, verbose = FALSE, ...)
    archs=aa$BY
    archslist=aa$run_results$archslist
  }
  #
  if(verbose){cat('Align archetypes by given final solution','\n')}
  #
  zz=align_archetypes_from_list(archslist,given_arch=archs,
                                varnames = colnames(df), verbose = verbose)
  dalist2=zz$archs_aligned
  names(dalist2)=paste0("after_iters_",1:length(dalist2))
  # Main computations
  if(doparallel){
    tp1=Sys.time()
    nworkers=nworkers
    #Function:
    runcloser=function(j,df,npoints,dalist2){
      archs=dalist2[[j]]
      jclosern=list()
      for(k in 1:dim(archs)[1]){
        mj=t(t(df)-c(unlist(archs[k,])))
        mj2=rowSums(mj^2)
        names(mj2)=1:length(mj2)
        jj=as.integer(names(sort(mj2)[1:npoints]))
        jclosern[[k]]=sort(jj)
      }
      closern=do.call(c,jclosern)
      return(closern)
    }
    environment(runcloser) <- .GlobalEnv
    cl <- makeCluster(nworkers);registerDoParallel(cl);
    clusterEvalQ(cl=cl,list(library("Matrix")));
    dclosern=parallel::parLapply(cl=cl,1:length(dalist2),runcloser,df=df,npoints=npoints,dalist2=dalist2)
    stopCluster(cl)
    tp2=Sys.time();
    if(verbose){cat(paste0('Time for computations in parallel was ',round(difftime(tp2,tp1,units="secs"),digits=3),' secs'),'\n')} 
  }else{
    #Function:
    runcloser=function(j,df,npoints,dalist2){
      archs=dalist2[[j]]
      jclosern=list()
      for(k in 1:dim(archs)[1]){
        mj=t(t(df)-c(unlist(archs[k,])))
        mj2=rowSums(mj^2)
        names(mj2)=1:length(mj2)
        jj=as.integer(names(sort(mj2)[1:npoints]))
        jclosern[[k]]=sort(jj)
      }
      closern=do.call(c,jclosern)
      return(closern)
    }
    t1=Sys.time()
    dclosern=lapply(1:length(dalist2),runcloser,df=df,npoints=npoints,dalist2=dalist2)
    t2=Sys.time()
    if(verbose){cat(paste0('Time for computations in serial was ',round(difftime(t2,t1,units="secs"),digits=3),' secs'),'\n')} 
  }
  #
  #
  dcloser=do.call(rbind,dclosern)
  mzn=data.frame(diff(as.matrix(dcloser)))
  mzn$sum=rowSums(mzn)
  #
  for(j in dim(mzn)[1]:2){
    zs=mzn$sum[j]
    if(zs!=0){
      jn=j+1
      if(verbose){print(jn)}
      break
    }
  } 
  #
  ptsfin=df[as.vector(tail(dcloser,1)),]
  #
  irdf=as.integer(rownames(ptsfin))
  #
  mpf=matrix(irdf,ncol=npoints,byrow=T)
  #
  out=list("rows_history"=dclosern,"iter_terminal"=jn,
           "rows_closer"=irdf,"rows_closer_matrix"=mpf,
           "solution_used"=aa)
  return(out)
}