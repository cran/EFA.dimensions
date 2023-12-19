
# parallel analysis of eigenvalues (no real data required)
PARALLEL <- function (Nvars=50, Ncases=300, Ndatasets=100, extraction='PCA', percentile=95,
                       corkind='pearson', verbose=TRUE, factormodel){
  
# deprecated  
if (!missing(factormodel))  extraction <- factormodel

evals <- matrix(0, nrow = Nvars, ncol = Ndatasets)
# pb <- utils::txtProgressBar(min = 0, max = Ndatasets, style = 3) 
for (nds in 1:Ndatasets) { 
	# utils::setTxtProgressBar(pb,nds)
#	message('Progress - Random Dataset: ', nds, 'of ', Ndatasets, '\r'); utils::flush.console()
	randat <- matrix(rnorm(Ncases*Nvars),nrow=Ncases,ncol=Nvars)
	# random data correlation matrix
	if (corkind=='pearson')      Rrand <- cor(randat, method='pearson') 
	if (corkind=='kendall')      Rrand <- cor(randat, method='kendall') 
	if (corkind=='spearman')     Rrand <- cor(randat, method='spearman') 
		
	# random data eigenvalues
	if (extraction=='PCA')  evals[,nds ] <- eigen(Rrand) $values 
	if (extraction=='PAF') {
		smc = 1 - (1 / diag(solve(Rrand)))
		diag(Rrand) <- smc
		evals[,nds] <- eigen(Rrand) $values 
	}
	if (extraction=='image') {
		d <-  diag(1 / diag(solve(Rrand)))
		gvv <- Rrand + d %*% solve(Rrand) %*% d - 2 * d
		s <- sqrt(d)                  #  Velicer 1974 p 565 formula (7)
		r2 <- solve(s) %*%  gvv  %*% solve(s)  # Velicer 1974 p 565 formula (5)
		evals[,nds] <- eigen(r2) $values
		} 	
} 
# mean & percentile eigenvalues for each position
means <- apply(evals, 1, mean) 
# sorting the eigenvalues for each root
for (luper in 1:Nvars) { evals[luper,] <- sort(evals[luper,]) }
percentiles <- as.matrix(evals[,round((percentile*Ndatasets)/100)])
results <- cbind(1:Nvars,means,percentiles)
rownames(results) <- matrix((''),nrow(results),1)
colnames(results) <- c('Root', 'Mean', 'Percentile')

if (verbose == TRUE) {
	message('\n\nPARALLEL ANALYSIS')
	# specification notices
	if (corkind=='pearson')    ctype <- 'pearson' 
	if (corkind=='kendall')    ctype <- 'kendall' 
	if (corkind=='spearman')   ctype <- 'spearman' 
	message('\nType of correlations specified for the random data eigenvalues: ', ctype)
	if (extraction=='PCA')    message('\nExtraction Method: Principal Components') 
	if (extraction=='PAF')    message('\nExtraction Method: Common Factor Analysis')
	if (extraction=='image')  message('\nExtraction Method: Image Factor Extraction')
	message('\nVariables = ', Nvars) 
	message('\nCases = ', Ncases) 
	message('\nNdatasets = ', Ndatasets) 
	message('\nPercentile = ', percentile, '\n') 
	print(round(results,3), print.gap=3)
}
parOutput <- list(eigenvalues=results)
return(invisible(parOutput))
}
