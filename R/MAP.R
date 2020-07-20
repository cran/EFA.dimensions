
MAP <- function (data, corkind='pearson', verbose=TRUE) {
# Velicer's MAP test -- takes raw data or a correlation matrix
data <- as.matrix(data)
     
nvars  <- ncol(data)
# determine whether data is a correlation matrix
if (nrow(data) == ncol(data)) {
	if (max(diag(data)) == 1 & min(diag(data)) == 1) {datakind = 'correlations'}} else{ datakind = 'notcorrels'}
if (datakind == 'correlations')  rdata <- data 
if (datakind == 'notcorrels') {
	ncases <- nrow(data)
	if (anyNA(data) == TRUE) {
		data <- na.omit(data)
		message('\nCases with missing values were found and removed from the data matrix.\n')
	}
	if (corkind=='pearson')     rdata <- cor(data, method='pearson') 
	if (corkind=='kendall')     rdata <- cor(data, method='kendall') 
	if (corkind=='spearman')    rdata <- cor(data, method='spearman') 
	if (corkind=='polychoric')  rdata <- POLYCHORIC_R(data) 
}
eigs    <- eigen(rdata)
eigval  <- diag(eigs$values)
eigvect <- eigs$vectors
# eigs  <- svd(rdata)
# eigval  <- diag(eigs$d)
# eigvect <- eigs$u
if (min(eigval) < 0) {
     message('\nThe correlation matrix is not positive definite, which means that')
     message('\nthere are eigenvalues less than zero. MAP test computations')
     message('\nare not possible.\n')
} else {
     
     rm(eigs)
     
     loadings <- eigvect %*% sqrt(eigval)
     
     fmfm4  <- matrix(NA,nvars,3)
     fmfm4[,1]  <- 0:(nvars-1)
     fmfm4[1,2] <- (sum(sum(rdata^2))-nvars)/(nvars*(nvars-1))
     fmfm4[1,3] <- (sum(sum(rdata^4))-nvars)/(nvars*(nvars-1))
     
     
     # pb <- utils::txtProgressBar(min = 0, max = (nvars - 1), style = 3) # create progress bar
     for (m in 1:(nvars - 1)) {
          #     Sys.sleep(0.1) # for the progress bar
          a <- loadings[,1:m]
          partcov <- as.matrix(rdata - tcrossprod(a,a))  # faster than as.matrix(rdata - (a %*% t(a)))
          
          if (max(partcov) > .0001) {
               d <- diag ( (1 / sqrt(diag(partcov))))
               pr <- d %*% (partcov %*% d)  # faster than d %*% partcov %*% d
               fmfm4[m+1,2] <- (sum(sum(pr^2))-nvars)/(nvars*(nvars-1))
               fmfm4[m+1,3] <- (sum(sum(pr^4))-nvars)/(nvars*(nvars-1))
          }	else {break}	
          
          #	rm(a,partcov,d,pr) # remove large matrices to free up memory re: R creates duplicates
          #     utils::setTxtProgressBar(pb, m) # update progress bar
     }
     # close(pb)
     
     
     # identifying the smallest fm values & their locations
     nfMAP   <- which.min(na.omit(fmfm4[,2])) - 1
     nfMAP4  <- which.min(na.omit(fmfm4[,3])) - 1
     
     dimnames(fmfm4) <-list(rep('', dim(fmfm4)[1]))
     colnames(fmfm4) <- c('root','  Avg.Corr.Sq.','  Avg.Corr.power4')
     
     evals <- cbind(1:nvars, (diag(eigval)))
     dimnames(evals) <-list(rep('', dim(evals)[1]))
     colnames(evals) <- c('root',' eigenvalue')
     
     if (verbose == TRUE) { 
          
          message("\n\nVelicer's Minimum Average Partial Test")
          
          if (datakind == 'correlations') message('\nThe entered data is a correlation matrix.') 
          
          if (datakind == 'notcorrels') {
               message('\nNumber of cases in the data file = ', ncases)
               message('\nNumber of variables in the data file = ', nvars)
               # message('\nSummary statistics for the data file:\n')
               # print(summary(data))
               
               # specification notices
               if (corkind=='pearson')    {message('\nCorrelations to be Analyzed: Pearson')}
               if (corkind=='kendall')    {message('\nCorrelations to be Analyzed: Kendall')}
               if (corkind=='spearman')   {message('\nCorrelations to be Analyzed: Spearman')}
               if (corkind=='polychoric') {message('\nCorrelations to be Analyzed: Polychoric')}
          }
          
          message('\nEigenvalues:\n')
          print(round(evals,5))
          
          message("\nVelicer's Average Squared Correlations\n")
          print(round(fmfm4,5))
          
          message('\nThe smallest average squared correlation is ', round(min(na.omit(fmfm4[,2])),5))
          message('\nThe smallest average 4rth power correlation is ', round(min(na.omit(fmfm4[,3])),5))
          
          message('\nThe Number of Factors According to the Original (1976) MAP Test is = ', nfMAP,  labels = NULL)
          message('\nThe Number of Factors According to the Revised (2000) MAP Test is = ', nfMAP4, labels = NULL, '\n')
     }
     
     mapOutput <- list(eigenvalues=evals, avgsqrs=fmfm4, nfMAP=nfMAP, nfMAP4=nfMAP4)
     
     return(invisible(mapOutput))
     
     message('\n')
}
}
