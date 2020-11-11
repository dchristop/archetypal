study_AAconvergence=function(df, kappas, method = "projected_convexhull", rseed = NULL,
                             chvertices = NULL, plot = TRUE, ...){
  # Run AA under default or chosen parameters (...)
  aa=archetypal(df = df, kappas = kappas, method = method, save_history = TRUE, rseed = rseed, ...)
  # Store results
  res=aa$run_results
  sse=res$SSE
  Blist=res$Blist
  archslist=res$archslist
  # Compute lowess of SSE
  ssel=lowess(1:length(sse),sse)
  # Find its UIK point
  uiklowess=ede(ssel$x,ssel$y,0)[1]
  ifinals=(uiklowess+1):length(sse)
  #############################################
  # Work for iterations after 'uiklowess' now:
  #############################################
  # Store final sequnece of B-matrices
  Bfinals=Blist[ifinals]
  names(Bfinals)=ifinals
  # Store final sequence of archetypes
  archsfinals=archslist[ifinals]
  names(archsfinals)=ifinals
  # Align final archetypes by using last one as guide
  archsaligned=align_archetypes_from_list(archsfinals,
                                          given_arch = archsfinals[[length(archsfinals)]],
                                          verbose = FALSE)$archs_aligned
  names(archsaligned)=ifinals
  #################################
  # Use Convex Hull vertices or not
  #################################
  ifelse(!is.null(chvertices),{ch=chvertices;hull=TRUE},{ch=NULL;hull=FALSE})
  # Percentage on CH, if asked
  ifelse(hull,{PCS=sapply(Bfinals, function(x){
    100*check_Bmatrix(x,chvertices = ch,verbose = FALSE)$used_on_convexhull})},{PCS=NA})
  ###############################################################################
  # Find Aitken extrapolations for lowess approximation of SSE's after 'uiklowess'
  ###############################################################################
  x=sse[ifinals]
  n=length(x)
  pxi=c()
  ei=c()
  for(i in 1:(n-2)){
    ei[i]=-(x[i+2]-x[i+1])^2/((x[i+2]-x[i+1])-(x[i+1]-x[i])) #Aitken  error formula
    pxi[i]=x[i+2]+ei[i] #Aitken extrapolation
  }
  ai=ifinals[1:(length(ifinals)-2)]
  daitken=data.frame("i"=ai,"pxi"=pxi,"ei"=ei)
  rownames(daitken)=ai
  ###############################
  # Order and Rate of convergence
  ###############################
  # Compute approximately the order of convergence p by using final SSE and the definitions of book:
  # Atkinson, K. E.,An Introduction to Numerical Analysis, Wiley & Sons,1989 : Definition of page 56 (62)
  # After taking natural logarithms and performing regression
  # Estimation of limit = last value of Aitken extrapolations
  a=tail(pxi,1)
  ################
  error_new=abs(a-x[2:n])
  error_old=abs(a-x[1:(n-1)])
  dlog=data.frame("log_error_new"=log(error_new),"log_error_old"=log(error_old))
  yy=lm(log_error_new~log_error_old, dlog)
  #
  #
  dcp=data.frame(coef(summary(yy)))
  colnames(dcp)=c("estimation","std.error","t.value","p.value")
  rownames(dcp)=c("log(c)","p")
  dcp
  #
  c_est=exp(dcp[1,1])
  p_est=dcp[2,1]
  #
  #
  ###############
  # Plot if asked
  ###############
  if(plot){
    #
    #
    def.par <- par(no.readonly = TRUE)
    par(mar = c(2.5, 3.5, 1.5, 0.5))
    par(mgp = c(1.5, 0.5, 0))
    mai = c(1, 1, 0.1, 0.1)
    par(oma = c(0, 0, 0, 0))
    #
    layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE),respect = FALSE)
    #1
    x1=1:length(sse)
    y1=sse
    plot(x1,y1,type='b',xlab="iteration",ylab="", pch=19,cex=0.5, 
         main = "UIK: all SSE", cex.axis = 0.8, tcl = -0.4, las = 1)
    knee1=uik(x1,y1)
    abline(v=knee1)
    text(x=knee1+2,y=max(y1), knee1, font = 2)
    legend("topright",pch=c(19,NA),lty=c(2,1),col=c('black','black'),
           legend=c("SSE","UIK"),bty="n")
    grid()
    #2
    plot(x1,y1,type='b',xlab="iteration", ylab="", pch = 19, cex = 0.5,
         main = "UIK: lowess(SSE)", cex.axis = 0.8, tcl = -0.4, las = 1)
    lines(ssel$x,ssel$y,col='blue',lwd=2,lty=1)
    abline(v=uiklowess,col='blue')
    text(x=uiklowess+2,y=max(y1),uiklowess, font = 2)
    legend("topright",pch=c(19,NA,NA),lty=c(4,1,1),lwd=c(1,2,1),
           col=c('black','blue','blue'),legend=c("SSE","lowess(SSE)","UIK"),bty="n")
    grid()
    #3
    x3=ifinals
    y3=sse[x3]
    plot(x3,y3,type='b',xlab="iteration", ylab="", pch=19,cex=0.5, 
         main="UIK: last SSE",cex.axis = 0.8, tcl = -0.4, las = 1)
    knee3=x3[ede(x3,y3,0)[1]]
    abline(v=knee3,col='black')
    text(x=knee3+2,y=max(y3),knee3, font = 2)
    legend("topright",pch=c(19,NA),lty=c(4,1),col=c('black','black'),
           legend=c("SSE","UIK"),bty="n")
    grid()
    #4
    if(hull){
      plot(daitken$i,daitken$pxi,type='b',pch=19,cex=0.5, xlab="iteration", ylab="",
           main="Aitken extrapolation SSE",cex.axis = 0.8, tcl = -0.4, las = 1)
      grid()
    }else{
      # extrapolation
      plot(daitken$i,daitken$pxi,type='b',pch=19,cex=0.5, xlab="iteration", ylab="",
           main="Aitken extrapolation SSE",cex.axis = 0.8, tcl = -0.4, las = 1)
      grid()
      # error
      plot(daitken$i,daitken$ei,type='b',pch=19,cex=0.5, xlab="iteration", ylab="",
           main="Aitken error SSE",cex.axis = 0.8,  tcl = -0.4, las = 1)
      grid()
    }
    #5
    if(hull){
      plot(ifinals,PCS,ylim=c(0,100),type='b',pch=19, xlab="iteration", ylab="",
           main="% of used rows on CH",cex.axis = 0.8,  tcl = -0.4, las = 1)
      grid()
      }
    #
    par(def.par)  
    #
    #
  }
  # Output
  out=list("SSE"=sse,"SSE_lowess"=ssel$y,"UIK_lowess"=uiklowess,"aitken"=daitken,
           "order_estimation"=p_est,"rate_estimation"=c_est,"significance_estimations"=dcp,
           "used_on_convexhull"=PCS, "aligned_archetypes"=archsaligned,
           "solution_used" = aa)
  ########
  # Return
  ########
  return(out)
}