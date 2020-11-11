align_archetypes_from_list=function (archs_list, given_arch = NULL, varnames = NULL, ndigits = 0, 
          parallel = FALSE, nworkers = NULL, verbose = TRUE) 
{
  ndims = prod(sapply(apply(sapply(archs_list, function(x) {
    dim(x)
  }), 1, unique), function(x) {
    length(x) == 1
  }))
  if (ndims == 0) {
    stop("The archetypes in archs_list must be of the same dimensions\n Please correct the list.")
  }
  if (is.null(varnames)) {
    varnames = colnames(archs_list[[1]])
  }
  ifelse(is.null(given_arch), {
    use_given_arch = FALSE
  }, {
    use_given_arch = TRUE
  })
  T1 = Sys.time()
  round2 <- function(x) {
    trunc(x + sign(x) * 0.5)
  }
  itrials = length(archs_list)
  dak = data.frame(do.call(rbind, archs_list))
  colnames(dak) = varnames
  k = dim(dak)[1]/itrials
  if (use_given_arch) {
    mp2 = NA
    arch_guide = given_arch
  } else {
    awords = apply(dak, 1, function(x, ndigits) {
      paste0(round2(x * (10^ndigits)), collapse = ":")
    }, ndigits)
    dak$word = awords
    mt = data.frame(table(awords))
    mt2 = mt[order(-mt$Freq), ]
    i1 = seq(1, length(awords), k) - 1
    j1 = i1 + 1
    j2 = i1 + k
    dphrases = lapply(1:itrials, function(i, j1, j2, awords) {
      sphrase = paste0(awords[j1[i]:j2[i]], collapse = "_")
      strial = rownames(dak[j1[i]:j2[i], ])
      return(cbind(phrase = as.character(sphrase), trial = strial))
    }, j1 = j1, j2 = j2, awords = awords)
    dphrases = data.frame(do.call(rbind, dphrases))
    mp = data.frame(table(dphrases$phrase))
    mp2 = mp[order(-mp$Freq), ]
    if(length(unique(mp2$Freq))==1){
      if(verbose){message(paste0("All phrases have equal frequency..."))}
      ii = sample(1:dim(mp2)[1],1)
      if(verbose){message(paste0("Lets choose in random the ... ",ii))}
    }else{
      mpmost = as.character(mp2[1, 1])
      imatch = which(dphrases$phrase %in% mpmost)
      mimatch=matrix(imatch,ncol=length(imatch)/k)
      #
      if(verbose){message(paste0("Most frequent phrases found at trials ... ",imatch))}
      ii = sample(mimatch[1,], 1)
      if(verbose){message(paste0("We choose most frequent which is ... ",ii))}
    }
    dakrows=ii:(ii+k-1)
    arch_guide = dak[dakrows, varnames]
  }
  if(verbose){
    message(paste0("Guide archetype will be next:"))
    print(arch_guide)
  }
  archguide = t(arch_guide)
  if (!parallel) {
    dlist = lapply(archs_list, function(archs, archguide) {
      sapply(as.list(data.frame(t(archs))), function(x, 
                                                     archguide) {
        sqrt(colSums((archguide - x)^2))
      }, archguide)
    }, archguide)
    align_list = lapply(dlist, function(x) {
      sol = lpSolve::lp.assign(x, "min")
      dd = sol$objval
      apply(sol$solution, 1, which.max)
    })
    archs_aligned = lapply(1:length(archs_list), function(i, 
                                                          archs_list, align_list) {
      data.frame(archs_list[[i]][align_list[[i]], ])
    }, archs_list, align_list)
  }
  else {
    if (is.null(nworkers)) {
      nc = detectCores()
      if (nc == 2) {
        nworkers = nc
      }
      else {
        nworkers = nc - 2
      }
    }
    tp1 = Sys.time()
    if(verbose){message(tp1)}
    fall = function(j, archs_list, archguide) {
      archs = archs_list[[j]]
      dj = sapply(as.list(data.frame(t(archs))), function(x, 
                                                          archguide) {
        sqrt(colSums((archguide - x)^2))
      }, archguide)
      sol = lpSolve::lp.assign(dj, "min")
      jalign = apply(sol$solution, 1, which.max)
      data.frame(archs[jalign, ])
    }
    environment(fall) <- .GlobalEnv
    cl <- makeCluster(nworkers)
    registerDoParallel(cl)
    clusterEvalQ(cl = cl, list(library("lpSolve")))
    archs_aligned = parallel::parLapply(cl = cl, 1:length(archs_list), 
                                        fall, archs_list, archguide)
    stopCluster(cl)
    tp2 = Sys.time()
    if(verbose){print(tp2 - tp1)}
  }
  outlist = list(arch_guide = arch_guide, phrases_most = mp2, 
                 archs_aa_output = dak, archs_aligned = archs_aligned)
  T2 = Sys.time()
  TT=T2-T1
  if(verbose){message(TT)}
  return(outlist)
}