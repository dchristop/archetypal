find_outmost_points=function(df,kappas){
  # Flatten data frame
  rownames(df)=1:dim(df)[1]
  # Find distances between rows:
  dm=as.matrix(dist(df,method = "euclidean", diag = F, upper = F, p = 2))
  # Find maximum distance per row and frequency table
  outmostrows=max.col(dm, ties.method="first")
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
    n1=kappas-nrowsmaxall
    r1=setdiff(1:dim(df)[1],rowsmaxall)
    rowsmaxall=c(rowsmaxall, sample(setdiff(1:dim(df)[1],rowsmaxall),n1))
  }
  #
  ########################################################################
  #Find the kappas outmost of outmost if possible
  rowsmaxkappas=as.integer(rowsmaxall[1:kappas])
  ########################################################################
  # Return
  out=list("outmost"=rowsmaxkappas,"outmostall"=rowsmaxall,"outmostfrequency"=di)
  return(out)
}