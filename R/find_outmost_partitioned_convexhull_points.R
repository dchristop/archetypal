find_outmost_partitioned_convexhull_points=function(df,kappas,np=10,nworkers=NULL){
  # Define number of workers if necessary
  if(is.null(nworkers)){
    nwall=parallel::detectCores()
    if(nwall<=2){nworkers=2}else{nworkers=nwall-2}
  }
  # Flatten data frame
  nn=dim(df)[1]
  rownames(df)=1:nn
  # Divide in approximately 'np' partitions:
  nk=round(nn/np) #equal sizes
  nr=nn-nk*(np-1) #rest size
  # Create np samples
  sizelist=list()
  irowsall=1:nn
  irows=c()
  nlist=c(rep(nk,np-1),nr)
  for( i in 1:np){    
    sizesample=sample(setdiff(irowsall,irows),nlist[i])
    irows=c(irows,sizesample)
    sizelist[[i]]=sizesample
  }
  cat(' ','\n')
  # Work for all sample in parallel mode:
  #
  tp1=Sys.time()
  #Function:
  runsizes=function(j,df,sizelist){
    # Find CH for each sample
    dhj=df[sizelist[[j]],]
    ch=as.list(convhulln(dhj,'Fx'))
    chj=unique(do.call(c,ch))
    chvertj=as.integer(rownames(dhj[chj,]))
    return(chvertj)
  }
  environment(runsizes) <- .GlobalEnv
  cl <- makeCluster(nworkers);registerDoParallel(cl);
  clusterEvalQ(cl=cl,list(library("geometry")));
  chsizes=parallel::parLapply(cl=cl,1:length(sizelist),runsizes,df=df,sizelist=sizelist)
  stopCluster(cl)
  tp2=Sys.time();print(tp2-tp1) 
  #
  # Create the overall data frame of sample champions
  #
  challsizes=do.call(c,chsizes)
  df2=df[challsizes,]
  print(dim(df2))
  nch=dim(df2)[1]
  pc=round(100*(nn-nch)/nn,2)
  # Work with it now
  if(pc<10){
    message(paste0('The reduction from initial data frame to ConvexHull data frame is ',pc,'% .
                 \n Since it is relatively small, the convex hull of convex hulls will not be computed.'))
    # Keep it as is
    chvertall=as.integer(rownames(df2))
    dfch=df2
    out=list("vertices"= chvertall,"points"= dfch)
    return(out)
  }else{
    # Find overall ConvexHUll
    ch=as.list(convhulln(df2,'Fx'))
    ch2=unique(do.call(c,ch))
    chvertall=as.integer(rownames(df2[ch2,]))
  }
  ########################################################################
  #Find data frame of CH rows
  rowschall=chvertall
  dfichall=df[chvertall,]
  nchall=length(rowschall)
  if(nchall==kappas){
    # We have found exactly the convex hull, no need to do anything else.
    outmostrows=as.integer(rownames(dfichall))
    di=data.frame("outmostrows"=outmostrows,"Freq"=rep(1,kappas),
                  "FreqPerCent"=rep(1/kappas,kappas),"CumFreqPerCent"=seq(1/kappas,1,by=1/kappas))
    rowsmaxkappas=outmostrows
    rowsmaxall=outmostrows
    out=list("outmost"=rowsmaxkappas,"outmostall"=rowsmaxall,"outmostfrequency"=di)
    return(out)
  }else{
    ########################################################################
    # Find outmost kappas
    #######################
    #
    # Find distances between rows:
    dm=as.matrix(dist(dfichall,method = "euclidean", diag = F, upper = F, p = 2))
    # Find maximum distance per row and frequency table
    outmostrows=max.col(dm, ties.method="first")
    outmostrows=as.integer(rownames(dfichall[outmostrows,]))
    di=as.data.frame(table(outmostrows),stringsAsFactors = F)
    di=di[order(-di$Freq),]
    di$outmostrows=as.integer(di$outmostrows)
    di$FreqPerCent=di$Freq/sum(di$Freq)
    di$CumFreqPerCent=cumsum(di$FreqPerCent)
    rownames(di)=1:dim(di)[1]
    rowsmaxall=di$outmostrows
    nrowsmaxall=length(rowsmaxall)
    #
    ####################################
    # Check if rows are less than kappas
    # Find the kappas outmost of outmost
    ####################################
    #
    #
    if(nrowsmaxall<kappas){
      if(nchall>=kappas){
        rowsmaxkappas=c(rowsmaxall,sample(setdiff(rowschall,rowsmaxall),kappas-nrowsmaxall))
      }else {
        n1=kappas-nchall
        r1=c(rowsmaxall,setdiff(rowschall,rowsmaxall))
        rowsmaxkappas=c(r1, sample(setdiff(1:dim(df)[1],r1),n1))
      }
    }else{
      rowsmaxkappas=rowsmaxall[1:kappas]
    }
    #
  }
  # Return
  out=list("outmost"=rowsmaxkappas,"outmostall"=rowsmaxall,"outmostfrequency"=di)
  return(out)
}