find_outmost_projected_convexhull_points=function(df,kappas,n=2){
  # Flatten data frame
  rownames(df)=1:dim(df)[1]
  #
  # Check given n
  #
  if(n>dim(df)[2]){
    warning(paste0('Given dimension of Projected Space was n = ',n,' was greater than column dimension of data frame\n It was set to n = ',dim(df)[2]))
    n=dim(df)[2]
  }
  #
  if(n==1){
    warning(paste0('Given dimension of Projected Space was n = ',n,' \n It was set to n = 2'))
    n=2
  }
  #
  # Find all possible pairs of variables-columns
  jcombs=combn(1:dim(df)[2],n)
  # Make list of convex hull for relevant facets
  if(n==dim(df)[2]){
    ch=as.list(convhulln(df,'Fx'))
    ichall=unique(do.call(c,ch))
  }else{
    ichlist=list()
    for(j in 1:dim(jcombs)[2]){
      if(n==2){
        # Use 'chull' for n=2 only
        ichlist[[j]]=chull(df[,jcombs[,j]])
      }else{
        chj=as.list(convhulln(df[,jcombs[,j]],'Fx'))
        ichlist[[j]]=unique(do.call(c,chj))
      }
    }
    # Store all of them
    ichall=unlist(ichlist)
  }
  # Find relevant vectors
  dfichall=data.frame(df[ichall,])
  # Store all available
  rowschall=unique(as.integer(rownames(dfichall)))
  nchall=length(rowschall)
  # Find distances between rows:
  dm=as.matrix(dist(dfichall,method = "euclidean", diag = F, upper = F, p = 2))
  # Find maximum distance per row abd frequency table
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
  ####################################
  #
  if(nrowsmaxall<kappas){
    if(nchall>kappas){
      rowsmaxall=c(rowsmaxall,sample(setdiff(rowschall,rowsmaxall),kappas-nrowsmaxall))
    }else{
      n1=kappas-nchall
      r1=c(rowsmaxall,setdiff(rowschall,rowsmaxall))
      rowsmaxall=c(r1, sample(setdiff(1:dim(df)[1],r1),n1))
    }
  }
  #
  ########################################################################
  #Find the kappas outmost of outmost
  rowsmaxkappas=as.integer(rowsmaxall[1:kappas])
  ########################################################################
  # Return
  out=list("outmost"=rowsmaxkappas,"outmostall"=rowsmaxall,"outmostfrequency"=di)
  return(out)
}
