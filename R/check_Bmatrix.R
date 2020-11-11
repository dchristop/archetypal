check_Bmatrix=function(B,chvertices=NULL,verbose=TRUE){
  #
  pw=list()
  ww=list()
  leadingr=c()
  leadingw=c()
  for(i in 1:dim(B)[1]){
    if(verbose){
      cat(paste0("Archetype ",i," is a mixture of only next rows with weights:"),"\n")
      cat(" ","\n")
    }
    dbi=B[i,]
    names(dbi)=1:dim(B)[2]
    rni=dbi[order(dbi,decreasing = T)]
    rni=rni[rni!=0]
    rowsi=as.integer(names(rni))
    leadingr[i]=rowsi[1]
    leadingw[i]=rni[1]
    pw[[i]]=rowsi
    ww[[i]]=rni
    if(verbose){
      print(rni)
      cat(" ","\n")
      cat("Weights add to:","\n")
      print(sum(rni))
      cat(" ","\n")
      cat(" ","\n")
    }
  }
  pwa=unlist(pw)
  if(!is.null(chvertices)){
    check=pwa%in%sort(chvertices)
    sumcheck=sum(pwa%in%sort(chvertices))
    pc=sumcheck/length(pwa)
    notok=!check
      if(verbose){
        cat('Used points that lie on convex hull are:','\n')
        cat(pwa[check],'\n')
        cat('Used points that do not lie on convex hull are:','\n')
        cat(pwa[notok],'\n')
        cat('Percentage of used points that really belong to convex hull is:','\n')
        cat(paste0(round(100*pc,2),' %'),'\n')
        cat(" ","\n")
      }
    }else{
      pc=NA
    }
 #
 #
 return(list("used_rows"=pwa,"used_weights"=ww,
             "leading_rows"=leadingr, "leading_weights"=leadingw,"used_on_convexhull"=pc))
}