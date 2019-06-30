archetypal=function(df,kappas,initialrows=NULL,method='projected_convexhull',nprojected=2,npartition=10,nfurthest=10,
                    maxiter=2000,conv_crit=1E-06,var_crit=0.9999,verbose=TRUE,rseed=NULL,
                    aupdate1=25,aupdate2=10,bupdate=10,muAup=1.2,muAdown=0.5,muBup=1.2,muBdown=0.5,
                    SSE_A_conv=1e-9,SSE_B_conv=1e-9,
                    save_history =FALSE,nworkers=NULL){
  # External Package usage: Matrix
  # Function that computes the PCHA for a data frame. 
  # It provides full control to the entire set of used parameters.
  # Based on Morten Morup's code http://www.mortenmorup.dk/index_files/Page327.htm , last accessed 2019-06-07
  #
  ################################
  # Internal update functions are:
  ################################
  #
  Aupdate=function(A,YYtBt,BYYtBt,muA,SST,SSE,niter=1,muAup=1.2,muAdown=0.5,nAup=1,nAdown=1,SSE_A_conv=1e-9){
    #Package usage: Matrix
    kappas=dim(A)[2]
    NN=dim(A)[1]
    ev=rbind(rep(1,kappas))
    #Iterations asked from input argument 'niter'
    for(k in 1:niter){
      # print(k)
      SSE_old=SSE
      gv=(A%*%BYYtBt-YYtBt)/(SST/NN)
      gv=gv-rowSums(A*gv)%*%ev
      stop=FALSE
      Aold=A
      #Main iterations
      itersmuDOWN=0 #index for preventing mu_Down to become zero
      while(!stop){
        A=Aold-gv*muA
        A[A<0]=0
        A=A/rowSums(A)
        # AtA=t(A)%*%A
        AtA=crossprod(A) #more efficient way for computing A'A
        SSE=SST-2*sum(A*YYtBt)+sum(AtA*BYYtBt)      
        #
        if(SSE<=SSE_old*(1+SSE_A_conv)){
          muA=muA*muAup
          stop=TRUE
          nAup=nAup+1 #counter for increasing mu
        }else{        
          muA=muA*muAdown
          nAdown=nAdown+1 #counter for decreasing mu
          itersmuDOWN=itersmuDOWN+1;if(itersmuDOWN>niter){stop=TRUE}
        }
        #
      }
    }
    #Return list of results
    return(list("A"=A,"SSE"=SSE,"muA"=muA,"AtA"=AtA,"nAup"=nAup,"nAdown"=nAdown))
  }
  #
  Bupdate=function(Y,AtY,BY,AtA,A,B,muB,SST,SSE,niter=1,muBup=1.2,muBdown=0.5,nBup=1,nBdown=1,SSE_B_conv=1e-9){
    #Package usage: Matrix
    eps=.Machine$double.eps
    N=dim(B)[2]
    et=rbind(rep(1,N))
    kappas=dim(B)[1]  
    AtYYt=tcrossprod(AtY, Y)  
    #Iterations asked from input argument 'niter'
    for(k in 1:niter){
      SSE_old=SSE   
      here1=tcrossprod(AtA%*%BY,Y)
      gt=(here1-AtYYt)/SST    
      gt=gt-rowSums(B*gt)%*%et
      stop2=FALSE
      B_old=B
      #Main iterations
      #
      while(!stop2){
        B=B_old-muB*gt
        B[B<0]=0
        nB=rowSums(B)+eps
        B=sparseMatrix(i=1:kappas, j=1:kappas, x = 1/nB,dims=c(kappas,kappas))%*%B
        B=as(B,"CsparseMatrix")
        #
        BY=B%*%Y 
        BYYtBt=BY%*%t(BY)
        YYtBt=Y%*%t(BY)
        SSE=SST-2*sum(A*YYtBt)+sum(AtA*BYYtBt)
        #Criterion of termination
        if(SSE<=SSE_old*(1+SSE_B_conv)){
          muB=muB*muBup
          stop2=TRUE
          nBup=nBup+1 #counter for increasing mu
        }else{        
          muB=muB*muBdown
          nBdown=nBdown+1 #counter for decreasing mu
        }
        #
      }
      #
    }
    #Return list of results
    return(list("B"=B,"SSE"=SSE,"muB"=muB,"BYYtBt"=BYYtBt,"BY"=BY,"nBup"=nBup,"nBdown"=nBdown))
    #
  }
  #    
  #######################################
  # Check if proper parameters are given:
  #######################################
  #
  if(conv_crit>1E-6){
    cat(paste0('Given conv_crit = ',conv_crit,' is relatively large, it was set to 1E-06'),'\n')
    conv_crit=1E-06
  }
  #
  if(var_crit<0.99){
    cat(paste0('Given var_crit = ',var_crit,' is relatively small, it was set to 0.9999'),'\n')
    var_crit=0.9999
  }
  #
  if(aupdate1<10){
    cat(paste0('Given aupdate1 = ',aupdate1,' is relatively small, it was set to 10'),'\n')
    aupdate1=10
  }
  #
  if(aupdate2<5){
    cat(paste0('Given aupdate1 = ',aupdate2,' is relatively small, it was set to 5'),'\n')
    aupdate2=5
  }
  #
  if(bupdate<5){
    cat(paste0('Given bupdate = ',bupdate,' is relatively small, it was set to 5'),'\n')
    bupdate=5
  }
  #
  if(muAup<1){
    cat(paste0('Given muAup = ',muAup,' cannot be less than 1, it was set to 1.1'),'\n')
    muAup=1.1
  }
  #
  if(muAup>2.5){
    cat(paste0('Given muAup = ',muAup,' is relatively large, it was set to 2.5'),'\n')
    muAup=2.5
  }
  #
  if(muBup<1){
    cat(paste0('Given muBup = ',muBup,' cannot be less than 1, it was set to 1.1'),'\n')
    muBup=1.1
  }
  #
  if(muBup>2.5){
    cat(paste0('Given muBup = ',muBup,' is relatively large, it was set to 2.5'),'\n')
    muBup=2.5
  }
  #
  if(muAdown > 0.5){
    cat(paste0('Given muAdown = ',muAdown,' is relatively large, it was set to 0.5'),'\n')
    muAdown =0.5
  }
  #
  if(muAdown < 0.1){
    cat(paste0('Given muAdown = ',muAdown,' is relatively small, it was set to 0.1'),'\n')
    muAdown=0.1
  }
  #
  if(muBdown > 0.5){
    cat(paste0('Given muBdown = ',muBdown,' is relatively large, it was set to 0.5'),'\n')
    muBdown =0.5
  }
  #
  if(muBdown < 0.1){
    cat(paste0('Given muBdown = ',muBdown,' is relatively small, it was set to 0.1'),'\n')
    muBdown=0.1
  }
  #
  if(SSE_A_conv>1E-9){
    cat(paste0('Given SSE_A_conv = ',SSE_A_conv,' is relatively large, it was set to 1E-09'),'\n')
    SSE_A_conv=1E-09
  }
  #
  if(SSE_B_conv>1E-9){
    cat(paste0('Given SSE_B_conv = ',SSE_B_conv,' is relatively large, it was set to 1E-09'),'\n')
    SSE_B_conv=1E-09
  }
  #
  ########################################
  # Define number of workers if necessary
  ########################################
  #
  if(is.null(nworkers)){
      nwall=parallel::detectCores()
      if(nwall<=2){nworkers=2}else{nworkers=nwall-2}
  }
  #
  ###################################################
  #Check if method='convexhull' is worth applying...
  ###################################################
  #
  if(method=='convexhull'){
  # Check if computing ConvexHull is time feasible and worth doing...
  vd=dim(df)[2]
  if(vd>6){
    cat(paste0("Since column dimension is ",vd," > 6, quick hull algorithm is not going to work fast. \n Now we change to method ='projected_convexhull' \n which will approximately perform the same task."),"\n")
    method='projected_convexhull'
    if(is.null(nprojected)){
      if(vd%in%2:3){
        nprojected=vd
      } else if(vd%in%4:5){
        nprojected=vd-1
      } else{
        nprojected=5
      }
    }
  }
  }
  #
  ##################
  # Begin main work
  ##################
  #
  T1=Sys.time() #overall time measure
  Y=as.matrix(df)
  NN=1:dim(Y)[1]
  #
  ########################################
  # Check if initial rows are given or not
  ########################################
  #
  if(is.null(initialrows)){
    #############################################################
    # Check if kappas=1 first: then use the point closest to mean 
    # Except method='furthestsum' where first row is chosen
    #############################################################
    if(kappas==1 & method!='furthestsum'){
      mdf=colMeans(df)
      dd=sqrt(rowSums((df-t(mdf))^2))
      irows=which.min(dd)
      freqstable=data.frame("outmostrows"=irows,"Freq"=1, "FreqPerCent"=1,"CumFreqPerCent"=1)
      YS=data.frame(df[irows,])
      cat('Next initial solution will be used...','\n')
      print(YS)
      cat(" ",'\n')
    }else{
      #
      ####################################
      # Find a good initial approximation
      ####################################
      #
      # Apply the given 'method' algorithm for initial solution finding
      #
      if(method=="projected_convexhull"){
        # Set maximum efficient n
        if(is.null(nprojected)){
          vd=dim(df)[2]
          if(vd%in%2:3){
            nprojected=vd
          } else if(vd%in%4:5){
            nprojected=vd-1
          } else{
            nprojected=5
          }
        }
        t1=Sys.time()
        projch=find_outmost_projected_convexhull_points(df,kappas = kappas,n=nprojected)
        t2=Sys.time();t12=round(as.numeric(t2-t1,units="secs"),digits=2) 
        message(paste0('Time for computing Projected Convex Hull was ',t12," secs"))
        irows=projch$outmost
        freqstable=projch$outmostfrequency
        #
        YS=Y[irows,]
        rownames(YS)=irows
        if(verbose){
          cat('Next projected convex hull initial solution will be used...','\n')
          print(YS)
          cat(" ",'\n')
        }
      } else if(method =="convexhull"){
        t1=Sys.time()
        ch=find_outmost_convexhull_points(df,kappas)
        t2=Sys.time();t12=round(as.numeric(t2-t1,units="secs"),digits=2) 
        message(paste0('Time for computing Convex Hull was ',t12," secs"))
        irows=ch$outmost
        freqstable=ch$outmostfrequency
        #
        YS=Y[irows,]
        rownames(YS)=irows
        if(verbose){
          cat('Next convex hull initial solution will be used...','\n')
          print(YS)
          cat(" ",'\n')
        }
      } else if(method =="partitioned_convexhull"){
        t1=Sys.time()
        parch=find_outmost_partitioned_convexhull_points(df,kappas,np=npartition,nworkers = nworkers)
        t2=Sys.time();t12=round(as.numeric(t2-t1,units="secs"),digits=2) 
        message(paste0('Time for computing Partitioned Convex Hull was ',t12," secs"))
        irows=parch$outmost
        freqstable=parch$outmostfrequency
        #
        YS=Y[irows,]
        rownames(YS)=irows
        if(verbose){
          cat('Next partitioned convex hull initial solution will be used...','\n')
          print(YS)
          cat(" ",'\n')
        }
      } else if(method =='furthestsum'){  
        # Check nfurthest
        if(nfurthest<10){
          message(paste0('Given nfurthest = ',nfurthest,' is relatively small, it was set to 10'))
          nfurthest=10
        }
        #
        t1=Sys.time() 
        fs=find_furthestsum_points(df=df,kappas=kappas,nfurthest=nfurthest,nworkers=nworkers,sortrows = TRUE)
        t2=Sys.time();t12=round(as.numeric(t2-t1,units="secs"),digits=2) 
        message(paste0('Time for computing ',nfurthest,' Furthest Sum initial solutions was ',t12," secs"))        
        irows=fs$outmost
        freqstable=fs$outmostfrequency
        #
        if(kappas!=1){YS=Y[irows,];rownames(YS)=irows}
        if(kappas==1){YS=data.frame(df[irows,])}
        if(verbose){
          cat('Next furthest fum initial solution will be used...','\n')
          print(YS)
          cat(" ",'\n')
        }
      }else if(method =="outmost"){
        # Find the kappas outmost points 
        t1=Sys.time() 
        yy=find_outmost_points(df=df,kappas = kappas)
        t2=Sys.time();t12=round(as.numeric(t2-t1,units="secs"),digits=2) 
        message(paste0('Time for computing outmost initial solutions for ',dim(df)[1],' rows was ',t12," secs"))        
        irows=yy$outmost
        freqstable=yy$outmostfrequency
        #
        YS=Y[irows,]
        rownames(YS)=irows
        if(verbose){
          cat('Next outmost initial solution will be used...','\n')
          print(YS)
          cat(" ",'\n')
        }
      } else if(method =="random"){
        #
        #############################################################################
        # If no method is given then use a random set of vectors as initial solution
        #############################################################################
        #
        irows=sample(NN,kappas)
        freqstable=NULL
        #
        YS=Y[irows,]
        rownames(YS)=irows
        if(verbose){
          cat('Next random initial solution will be used...','\n')
          print(YS)
          cat(" ",'\n')
        }
      } else {
        stop("You must give one of the following 'method' names: \n 'projected_convexhull', 'convexhull', 'partitioned_convexhull',
             'furthestsum', 'outmost', 'random'")
      }
      #
    }
  } else{ 
    #
    #############################
    # Use the given initial rows
    #############################
    #
    irows=initialrows
    freqstable=data.frame("outmostrows"=initialrows,"Freq"=1, "FreqPerCent"=1,"CumFreqPerCent"=1)
    #
    YS=Y[irows,]
    rownames(YS)=irows
    if(verbose){
      print(irows)
      cat('The initial solution that will be used is','\n')
      print(YS)
      cat(" ",'\n')
    }
  }
  # 
  ###############################
  # Proceed to main PCHA code
  ###############################
  #
  t1=Sys.time() #measure time for initial 'aupdate1' A updates
  #Compute SST
  SST=sum(Y^2)
  #Define C
  irows2=1:length(irows)
  ij=cbind(irows,irows2)
  Bt=sparseMatrix(i=ij[,1], j=ij[,2], x = rep(1,kappas),dims=c(length(NN),kappas))
  B=t(Bt)
  #Create initial BY from FurthestSum:
  BY=B%*%Y
  # Store it
  initialsolution=as.matrix(BY)
  rownames(initialsolution)=irows
  #
  ################
  #
  # Initialize mus:  
  muA=1
  muB=1
  # counters of increasing or decreasing mu
  nAup=1;nAdown=1
  nBup=1;nBdown=1
  # Begin loop
  # Initialize A 
  YYtBt=Y%*%t(BY)
  BYYtBt=BY%*%t(BY)
  #
  #Define and use seed for A
  if(!is.null(rseed)){
    set.seed(rseed);rnumbers=runif(kappas*length(NN))
    A=-log(matrix(rnumbers,length(NN),kappas,byrow = T))
  }else{
    A=-log(matrix(runif(kappas*length(NN)),length(NN),kappas,byrow = T))
  }
  A=A/rowSums(A)
  # AtA=t(A)%*%A
  AtA=crossprod(A) #efficient way of computing
  SSE=SST-2*sum(A*YYtBt)+sum(AtA*BYYtBt)
  ev=rbind(rep(1,kappas))
  #
  # Initial 'aupdate1' A updates
  #
  ya1=Aupdate(A,YYtBt,BYYtBt,muA,SST,SSE,niter=aupdate1,muAup=muAup,muAdown=muAdown,nAup=nAup,nAdown=nAdown,SSE_A_conv=SSE_A_conv)
  #Update values
  A=ya1$A
  SSE=ya1$SSE
  muA=ya1$muA
  AtA=ya1$AtA
  nAup=ya1$nAup
  nAdown=ya1$nAdown
  # 
  t2=Sys.time();t12=round(as.numeric(t2-t1,units="secs"),digits=2) #print time for initial 'aupdate1' A updates
  if(verbose){cat(paste0("Time for the ",aupdate1," initial A updates was ",t12," secs"),'\n')}
  # Define print layout
  dheader = sprintf('%12s | %12s | %15s | %10s | %10s | %10s | %10s | %12s | %12s','Iteration','VarExpl.','SSE','|dSSE|/SSE','muB','muA',' Time(s) ','nAup_nAdown','nBup_nBdown');
  dline = sprintf('|------------|--------------|-----------------|------------|------------|------------|------------|--------------|-------------|')
  #
  ###########################################
  # Set parameters for main iteration section
  ###########################################
  #
  iter=0
  dSSE=Inf
  varexplv=(SST-SSE)/SST
  #
  # Check if initial method has already found the optimal solution:
  #
  if(varexplv>var_crit){
    dSSE=0
    if(verbose){
      cat(paste0('Optimal solution was found from FurthestSum because Variance Explained = ',round(varexplv,6),' > ',var_crit),'\n')
      cat(" ",'\n')
      cat(dheader,"\n");cat(dline,"\n")
      cat(sprintf('%12.0f | %12.6f | %15.6e | %10.2e | %10.2e | %10.2e | %10.2f | %12s | %12s \n',iter,varexplv,SSE,abs(dSSE)/SSE,muB,muA,t1-t1,paste0(c(nAup,nAdown),collapse = "_"),paste0(c(nBup,nBdown),collapse = "_")))
      cat(dline,"\n")
    }
    # Sort components according to their importance
    csums=colSums(A)
    names(csums)=1:length(csums)
    isortv=as.integer(names(sort(csums,decreasing = T)))
    A=as.matrix(A[,isortv])
    B=as.matrix(B[isortv,])
    BY=as.matrix(BY[isortv,])
    if(verbose){
      cat(" ",'\n')
      cat(" BY = ",'\n')
      print(BY)
      cat(" ",'\n')
    }
    #
    # Set proper run results 
    #
    if(save_history ){
      run_results=list("SSE"=SSE,"varexpl"=varexplv,"time"=0,"Blist"=list(B),"archslist"=list(data.frame(as.matrix(BY))))
    }else{
      run_results=NULL
    }
    #
    T2=Sys.time();T12=round(as.numeric(T2-T1,units="secs"),digits=2) #print overall function execution time 
    # Return list of results
    return(list("BY"=BY,"A"=A,'B'=B,"SSE"=SSE,"varexpl"=varexplv,
                "initialsolution"=initialsolution,"freqstable"=freqstable,
                "iterations"=iter,"time"=T12,"converges"=TRUE,
                "nAup"=nAup,"nAdown"=nAdown,"nBup"=0,"nBdown"=0,
                "run_results"=run_results))
  }
  #
  ##################
  # Main iterations
  ##################
  #
  # Define run results if asked so, otherwise set it nulll
  #
  if(save_history){
    vsse=c()
    vvarexpl=c()
    vtime=c()
    Blist=list()
    archslist=list()
  }else{
    run_results=NULL
  }
  #
  #
  if(verbose){cat(dline,"\n");cat(dheader,"\n");cat(dline,"\n")}
  while(abs(dSSE)>=conv_crit*abs(SSE) & iter<maxiter & varexplv<var_crit){
    t1=Sys.time() #count iteration time
    iter=iter+1
    SSEold=SSE
    #
    # Store counters before the two B and A update rounds begin:
    #
    nBup_old=nBup;nBdown_old=nBdown
    nAup_old=nAup;nAdown_old=nAdown
    #
    # 
    ######################
    # AtY=t(A)%*%Y
    AtY=crossprod(A,Y) #efficient way of computing
    #
    # Update B 
    #
    yb=Bupdate(Y,AtY,BY,AtA,A,B,muB,SST,SSE,niter=bupdate,muBup=muBup,muBdown=muBdown,nBup=nBup,nBdown=nBdown,SSE_B_conv=SSE_B_conv)
    #
    # Update values
    #
    B=yb$B
    SSE=yb$SSE
    muB=yb$muB
    BYYtBt=yb$BYYtBt
    BY=yb$BY
    nBup=yb$nBup
    nBdown=yb$nBdown
    #
    # Update A
    #
    YYtBt=Y%*%t(BY)
    ya2=Aupdate(A,YYtBt,BYYtBt,muA,SST,SSE,niter=aupdate2,muAup=muAup,muAdown=muAdown,nAup=nAup,nAdown=nAdown,SSE_A_conv=SSE_A_conv)
    #
    # Update values
    #
    A=ya2$A
    SSE=ya2$SSE
    muA=ya2$muA
    AtA=ya2$AtA
    nAup=ya2$nAup
    nAdown=ya2$nAdown
    #
    # Evaluate dSSE and display iteration details if asked so
    # 
    dSSE=SSEold-SSE
    t2=Sys.time();t12=round(as.numeric(t2-t1,units="secs"),digits=2) #compute iteration time in secs
    #
    varexplv=(SST-SSE)/SST
    #
    # Store history of run results if asked
    #
    if(save_history ){
      vsse=c(vsse,SSE)
      vvarexpl=c(vvarexpl,varexplv)
      vtime=c(vtime,t12)
      Blist[[iter]]=B
      archslist[[iter]]=data.frame(as.matrix(BY))
    }
    #
    if(verbose){
      cat(sprintf('%12.0f | %12.9f | %15.9e | %10.2e | %10.2e | %10.2e | %10.2f | %12s | %12s \n',iter,varexplv,SSE,abs(dSSE)/SSE,muB,muA,t12,
                  paste0(c(nAup-nAup_old,nAdown-nAdown_old),collapse = "_"),paste0(c(nBup-nBup_old,nBdown-nBdown_old),collapse = "_")))
    }
    #
  }
  #
  if(verbose){cat(dline,"\n")}
  # Sort components according to their importance by inspecting weights of matrix A
  csums=colSums(A)
  names(csums)=1:length(csums)
  isort=as.integer(names(sort(csums,decreasing = T)))
  A=as.matrix(A[,isort])
  #
  if(kappas==1){BY=t(as.matrix(BY[isort,]));B=as.matrix(B)}else{BY=as.matrix(BY[isort,]);B=as.matrix(B[isort,])}
  #
  if(verbose){
    cat(" ",'\n')
    cat(" BY = ",'\n')
    print(BY)
    cat(" ",'\n')
  }
  #
  # Check if convergence exists:
  ##############################
  #
  if(iter<maxiter){iflag=TRUE}else{iflag=FALSE}
  #
  # Create list of run_results if asked so, otherwise leave it null
  if(save_history ){run_results=list("SSE"=vsse,"varexpl"=vvarexpl,"time"=vtime,"Blist"=Blist,"archslist"=archslist)}
  #
  ##########
  # Return:
  ##########
  #
  T2=Sys.time();T12=round(as.numeric(T2-T1,units="secs"),digits=2) #compute overall function execution time
  # Return list of results
  return(list("BY"=BY,"A"=A,'B'=B,"SSE"=SSE,"varexpl"=varexplv,
              "initialsolution"=initialsolution,"freqstable"=freqstable,
              "iterations"=iter,"time"=T12,"converges"=iflag,
              "nAup"=nAup,"nAdown"=nAdown,"nBup"=nBup,"nBdown"=nBdown,
              "run_results"=run_results))
  # 
}