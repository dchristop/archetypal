find_pcha_optimal_parameters<-function(df, kappas, method = "projected_convexhull", testing_iters = 10, 
                                       nworkers = NULL, nprojected = 2, npartition = 10, nfurthest = 100,
                                       sortrows = FALSE, mup1 = 1.1, mup2 = 2.50, mdown1 = 0.1, mdown2 = 0.5, 
                                       nmup = 10, nmdown = 10, rseed = NULL, plot = FALSE, ...){
  #######################################################################################
  # Create a 2D grid in parameter spaces (muA_up, muA_down) and (muB_up, muB_down)
  # Run archetypal on it and find 'ceteris paribus' the optimal pair of parameters
  # Criterion: smallest SSE found on last iteration
  #######################################################################################
  #
  ###################
  # Set seed if given
  ###################
  #
  if(!is.null(rseed)){rseed=rseed}else{rseed=NULL}
  #
  # Define number of workers if necessary
  #
  if(is.null(nworkers)){
    nwall=parallel::detectCores()
    if(nwall<=2){nworkers=nwall}else{nworkers=nwall-2}
  }
  #
  ######################################################
  # CREATE THE GRID OF PARAMETER SPACES (mu_up, mu_down)
  ######################################################
  #
  mu_ups=seq(mup1,mup2,len=nmup)
  mu_downs=seq(mdown1,mdown2,len=nmdown)
  mu_grid=expand.grid(mu_ups,mu_downs)
  colnames(mu_grid)=c("mu_up","mu_down")
  #
  ###############################################
  # FIND INITIAL APPROXIMATION GIVENT THE method
  ###############################################
  #
  #
  if(method=="projected_convexhull"){
    sol_initial=find_outmost_projected_convexhull_points(df = df, kappas = kappas, npr = nprojected)$outmost
  }else if(method=="convexhull"){
    sol_initial=find_outmost_convexhull_points(df = df, kappas = kappas)$outmost
  }else if(method=="partitioned_convexhull"){
    sol_initial=find_outmost_partitioned_convexhull_points(df = df, kappas = kappas, np = npartition, nworkers = nworkers)$outmost
  }else if(method=="furthestsum"){
    sol_initial=find_furthestsum_points(df =df, kappas = kappas, nfurthest = nfurthest, sortrows = sortrows)$outmost
  }else if(method=="outmost"){
    sol_initial=find_outmost_points(df = df, kappas = kappas)$outmost
  }else if(method=="random"){
    sol_initial=sample(1:dim(df)[1],kappas, replace = FALSE)
  }else{
    stop("Please use as method = ' ' one of next available methods: \n 
         'projected_convex_hull', 'convexhull', 'partitioned_convexhull', 'outmost', 'random'")
  }
  #
  #####################
  # WORK IN PARALLEL
  #####################
  # 
  run_progress=list()
  # Function:
  runaa=function(i,k,df,sol_initial,testing_iters,mu_grid,rseed,...){
    mu_up=mu_grid[i,1]
    mu_down=mu_grid[i,2]
    yy=archetypal(df = df, kappas = k, initialrows = sol_initial, maxiter = testing_iters,
                  muAup = mu_up, muAdown = mu_down, muBup = mu_up, muBdown = mu_down,
                  verbose = FALSE, rseed = rseed,
                  save_history = TRUE)
    list("SSE_iters"=yy$run_results$SSE,"varexpl_iters"=yy$run_results$varexpl,
         "TIME_iters"=yy$run_results$time,"SSE"=yy$SSE,"varexpl"=yy$varexpl)
  }
  environment(runaa) <- .GlobalEnv
  tp1=Sys.time()
  cl <- makeCluster(nworkers);registerDoParallel(cl);
  clusterEvalQ(cl=cl,list(library("archetypal")))
  run_progress=parallel::parLapply(cl=cl,1:dim(mu_grid)[1],runaa,
                                   k=kappas,df=df,sol_initial=sol_initial,
                                   testing_iters=testing_iters, mu_grid=mu_grid,
                                   rseed=rseed,...)
  stopCluster(cl)
  tp2=Sys.time();cat(" ","\n");print(tp2-tp1) 
  #
  ########################
  # Find minimum final SSE
  ########################
  #
  lsse=sapply(run_progress, function(x){x$SSE})
  imin=which.min(lsse)
  minsse=lsse[imin]
  #
  # Set optimal
  #
  mu_up_opt=mu_grid[imin,"mu_up"]
  mu_down_opt=mu_grid[imin,"mu_down"]
  outbasics=c(mu_up_opt,mu_down_opt,minsse)
  names(outbasics)=c("mu_up_opt","mu_down_opt","min_sse")
  print(outbasics)
  #
  #################
  # PLOT (if asked)
  #################
  #
  # Compute the linear regression (z = ax + by + d) or full quadratic regression
  # fit <- lm(z ~ x + y)
  x=mu_grid[,1];y=mu_grid[,2];z=lsse
  fit <- lm(z ~ (x+y)^2+I(x^2)+I(y^2))
  # predict values on regular xy grid
  grid.lines = 30
  x.pred <- seq(min(x), max(x), length.out = grid.lines)
  y.pred <- seq(min(y), max(y), length.out = grid.lines)
  xy <- expand.grid( x = x.pred, y = y.pred)
  z.pred <- matrix(predict(fit, newdata = xy), 
                   nrow = grid.lines, ncol = grid.lines)
  # fitted points for droplines to surface
  fitpoints <- predict(fit)
  # scatter plot with regression plane
  maintitle=paste0("mu_up = ",round(mu_up_opt,2)," , mu_down = ",round(mu_down_opt,2),
                   "\n iters = ",testing_iters," , SSE [x] = ",round(minsse,3))
  scatter3D(x, y, z,pch = 18, cex = 2,colkey=F,colkeyside = 1, 
            main = maintitle, cex.main=0.75,
            theta = 20, phi = 30, ticktype = "detailed", bty ="g",
            xlab = "mu_up", ylab = "mu_down", zlab = "SSE",  
            surf = list(x = x.pred, y = y.pred, z = z.pred,  
                        facets=NA, fit = fitpoints), colvar = z, cex.axis=0.8, cex.lab=0.7)
  scatter3D(x[imin],y[imin],z[imin],colkey=F,add=TRUE,col='black',cex=4,pch=7)
  #
  #
  out=list("mu_up_opt"=mu_up_opt,"mu_down_opt"=mu_down_opt,"min_sse"=minsse,"seed_used"=rseed,
           "sol_initial"=sol_initial,"method_used"=method,"testing_iters"=testing_iters)
  #
  ########
  # Return
  ########
  #
  return(out)
}