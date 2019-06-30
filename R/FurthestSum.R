FurthestSum=function(Y,kappas,irows,exclude=NULL){
  #Package usage: Matrix
  #Initialize
  N=dim(Y)[1]  
  index=1:N;
  index[exclude]=0;
  index[irows]=0
  irun=irows
  sum_distance=rep(0,N)
  #
  Y=as.matrix(Y)
  Y2=rowSums(Y^2)
  # Main loop
  for(k in 1:(kappas+10)){
    if(k>(kappas-1)){
      #Remove initial seed 
      Yq=as.vector(Y%*%Y[irows[1],])
      sum1=as.vector(Y2)-2*Yq+rep(Y2[irows[1]],length(Yq))
      sum1[sum1<0]=0
      sum_distance=sum_distance-sqrt(sum1) 
      index[irows[1]]=irows[1]
      if(length(irows)==1){irows=c()}else{irows=irows[2:length(irows)]}
    }
    inonzero=which(index!=0)
    Yq=as.vector(Y%*%matrix(Y[irun,],ncol=1))
    sum2=as.vector(Y2)-2*Yq+rep(Y2[irun],length(Yq))
    sum2[sum2<0]=0
    sum_distance=sum_distance+sqrt(sum2)
    indt=which.max(sum_distance[inonzero])
    if(length(indt)==0){return(irows)}
    irun=inonzero[indt]
    irows=c(irows,irun)
    index[irun]=0
  }
  return(irows)
}