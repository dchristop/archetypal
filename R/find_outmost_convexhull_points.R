find_outmost_convexhull_points=function(df, kappas){
  # Flatten data frame
  rownames(df)=1:dim(df)[1]
  # Find convex hull of data frame
  ch=as.list(convhulln(df,'Fx'))
  rowschall=unique(do.call(c,ch))
  dfichall=data.frame(df[rowschall,])
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