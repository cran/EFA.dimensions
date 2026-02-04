

RAWPAR <- function (data, randtype='generated', extraction='PCA', 
                    Ndatasets=100, percentile=95,
                    corkind='pearson', corkindRAND='pearson', 
                    Ncases=NULL, eval_plot=TRUE, verbose=TRUE, factormodel) {
  
  # deprecated  
  if (!missing(factormodel))  extraction <- factormodel
  
  data <- MISSING_DROP(data)
  
  Nvars  <- ncol(data)
  
  # set up cormat
  cordat <- setupcormat(data, corkind=corkind, Ncases=Ncases)
  R_real <- cordat$cormat
  ctype  <- cordat$ctype
  Ncases <- cordat$Ncases
  datakind <- cordat$datakind
  
  # determine whether the data are whole numbers
  if (all((data - round(data)) == 0)) { wholenums=1 } else { wholenums=0 }
  
  # real data eigenvalues
  realevals <- eigvals(cormatrix=R_real, extraction=extraction) 
  
  # random data
  if (datakind == 'correlations')  randtype = 'generated' 
  evals <- matrix(0, nrow = Nvars, ncol = Ndatasets)
  # pb <- utils::txtProgressBar(min = 0, max = Ndatasets, style = 3) 
  for (nds in 1:Ndatasets) {
    
    # utils::setTxtProgressBar(pb,nds)
    # message('Progress - Random Dataset: ', nds, 'of ', Ndatasets, '\r'); utils::flush.console()
    
    # random data -- generated or permuted 
    if (randtype == 'generated') randat <- matrix(rnorm(Ncases*Nvars),nrow=Ncases,ncol=Nvars) 
    if (randtype == 'permuted')  randat <- apply(data, 2, sample) 
    
    # random data correlation matrix
    if (corkindRAND=='pearson')     { R_rand <- cor(randat, method='pearson') }
    if (corkindRAND=='kendall')     { R_rand <- cor(randat, method='kendall') }
    if (corkindRAND=='spearman')    { R_rand <- cor(randat, method='spearman') }
    
    if (corkindRAND=='polychoric' & wholenums==1)  { R_rand <- POLYCHORIC_R(randat, verbose=FALSE) }
    if (corkindRAND=='polychoric' & wholenums==0)  { R_rand <- cor(randat, method='pearson') }
    if (corkindRAND=='gamma' & wholenums==1)       { R_rand <- Rgamma(randat, verbose=FALSE) }
    if (corkindRAND=='gamma' & wholenums==0)       { R_rand <- cor(randat, method='pearson') }
    
    # random data eigenvalues
    evals[,nds] <- eigvals(cormatrix=R_rand, extraction=extraction) 
  } 
  
  # mean & percentile eigenvalues for each position
  means <- apply(evals, 1, mean) 
  # sorting the eigenvalues for each root
  for (luper in 1:Nvars) { evals[luper,] <- sort(evals[luper,]) }
  percentiles <- as.matrix(evals[,round((percentile*Ndatasets)/100)])
  realevals <- as.matrix(realevals)
  NfactorsPA <- 0
  for (root in 1:Nvars) { 
    if (realevals[root,1] < percentiles[root,1]) { NfactorsPA <- root - 1; break } 
  }
  
  results <- cbind(1:Nvars,realevals,means,percentiles)
  rownames(results) <- matrix((''),nrow(results),1)
  colnames(results) <- c('Root', 'Real_Data', 'Mean', 'Percentile')
  
  results <- as.data.frame(results) 
  
  if (eval_plot) {
    
    plot(results$Root, results$Real_Data, pch=15, 
         xlab='Root', ylab='Eigenvalue', main='Eigenvalues Plot for Parallel Analyses',
         cex=0.8, cex.axis = 1, cex.lab=1, col='blue', type='b')
    axis(1, at=results$Root)
    
    lines(results$Root, results$Percentile, lwd=1, lend=1, col='red', cex=0.8, 
          type='b', lty=1, pch = 16)
    
    legend('topright', legend=c("Real Data", "PA percentile"),
           col=c("blue","red"), lty=1, cex=1, pch=c(15,16)) #, bty = "n")
    
    grid(nx = NULL, ny = NULL, lty = 2, col = "gray", lwd = .5)
  }
  
  
  if (verbose == TRUE) {
    
    message('\n\nPARALLEL ANALYSIS')
    if (datakind == 'correlations') message('\nThe entered data is a correlation matrix.') 
    if (randtype == 'generated')    message('\nRandomization method: generated data') 
    if (randtype == 'permuted')     message('\nRandomization method: permuted data') 
    
    message('\nType of correlations specified for the real data eigenvalues: ', ctype)
    
    message('\nType of correlations specified for the random data eigenvalues: ', corkindRAND)
    
    if (randtype == 'permuted' & datakind == 'correlations') {
      message('\nThe permuted data option does not work when the entered data are a correlation matrix.')
      message('\nSwitching to generated random data instead of permuted random data.')
    }
    
    if (extraction=='PCA')    message('\nExtraction Method: Principal Components') 
    if (extraction=='PAF')    message('\nExtraction Method: Common Factor Analysis')
    if (extraction=='image')  message('\nExtraction Method: Image Factor Extraction') 
    
    message('\nVariables = ', Nvars) 
    message('\nCases = ', Ncases) 
    message('\nNdatasets = ', Ndatasets) 
    message('\nPercentile = ', percentile) 
    message('\nCompare the Real Data eigenvalues below to the corresponding')
    message('random data Mean and Percentile eigenvalues\n\n')
    print(round(results,3), print.gap=3, row.names = FALSE)
    message('\nThe number of factors according to the parallel analysis = ', NfactorsPA, '\n')
  }
  
  rawparOutput <- list(eigenvalues=results, NfactorsPA=NfactorsPA)
  return(invisible(rawparOutput))
}
