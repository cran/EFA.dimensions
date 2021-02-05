
RAWPAR <- function (data, randtype='generated', factormodel='PCA', Ndatasets=100, percentile=95,
                    corkind='pearson', corkindRAND='pearson', Ncases=NULL, verbose=TRUE){

Nvars  <- ncol(data)

# set up cormat
cordat <- setupcormat(data, corkind=corkind, Ncases=Ncases)
cormat <- cordat$cormat
ctype  <- cordat$ctype
Ncases <- cordat$Ncases
datakind <- cordat$datakind


# determine whether the data are whole numbers
if (all((data-round(data)) == 0) == TRUE) { wholenums=1 
} else { wholenums=0 }
# real data correlation matrix
if (datakind == 'correlations') Rreal <- data
if (datakind == 'notcorrels')   { 
	Ncases <- nrow(data) 	
	if (corkind=='pearson')         { Rreal <- cor(data, method='pearson') }
	if (corkind=='kendall')         { Rreal <- cor(data, method='kendall') }
	if (corkind=='spearman')        { Rreal <- cor(data, method='spearman') }
	if (corkind=='polychoric')      { 
		if (wholenums == 1) { Rreal <- POLYCHORIC_R(data) } 
		if (wholenums == 0) { Rreal <- cor(data, method='pearson') }}
	if (corkindRAND=='gamma' & wholenums==1)       { Rreal <- Rgamma(randat, verbose=FALSE) }
	if (corkindRAND=='gamma' & wholenums==0)       { Rreal <- cor(randat, method='pearson') }
}

# real data eigenvalues
if (factormodel=='PCA')  realevals <- eigen(cormat) $values 
if (factormodel=='PAF') {
	smc = 1 - (1 / diag(solve(cormat)))
	diag(cormat) <- smc
	realevals <- eigen(cormat) $values 
}
if (factormodel=='image') { # Gorsuch 1983, p 113; Velicer 1974, EPM, 34, 564
	d <-  diag(1 / diag(solve(cormat)))
	gvv <- cormat + d %*% solve(cormat) %*% d - 2 * d
	s <- sqrt(d)                  #  Velicer 1974 p 565 formula (7)
	r2 <- solve(s) %*%  gvv  %*% solve(s)  # Velicer 1974 p 565 formula (5)
	realevals <- cbind(eigen(r2) $values) 
}

# random data
if (datakind == 'correlations')  randtype = 'generated' 
evals <- matrix(0, nrow = Nvars, ncol = Ndatasets)
# pb <- utils::txtProgressBar(min = 0, max = Ndatasets, style = 3) 
for (nds in 1:Ndatasets) { 
# 		utils::setTxtProgressBar(pb,nds)
	#message('Progress - Random Dataset: ', nds, 'of ', Ndatasets, '\r'); utils::flush.console()
	# random data -- generated or permuted 
	if (randtype == 'generated') randat <- matrix(rnorm(Ncases*Nvars),nrow=Ncases,ncol=Nvars) 
	if (randtype == 'permuted')  randat <- apply(data, 2, sample) 
	# random data correlation matrix
	if (corkindRAND=='pearson')     { Rrand <- cor(randat, method='pearson') }
	if (corkindRAND=='kendall')     { Rrand <- cor(randat, method='kendall') }
	if (corkindRAND=='spearman')    { Rrand <- cor(randat, method='spearman') }
	if (corkindRAND=='polychoric' & wholenums==1)  { Rrand <- POLYCHORIC_R(randat) }
	if (corkindRAND=='polychoric' & wholenums==0)  { Rrand <- cor(randat, method='pearson') }
	if (corkindRAND=='gamma' & wholenums==1)       { Rrand <- Rgamma(randat, verbose=FALSE) }
	if (corkindRAND=='gamma' & wholenums==0)       { Rrand <- cor(randat, method='pearson') }
		
	# random data eigenvalues
	if (factormodel=='PCA') { evals[,nds ] <- eigen(Rrand) $values }
	if (factormodel=='PAF') {
		smc = 1 - (1 / diag(solve(Rrand)))
		diag(Rrand) <- smc
		evals[,nds] <- eigen(Rrand) $values }
	if (factormodel=='image') {
		d <-  diag(1 / diag(solve(Rrand)))
		gvv <- Rrand + d %*% solve(Rrand) %*% d - 2 * d
		s <- sqrt(d)                  #  Velicer 1974 p 565 formula (7)
		r2 <- solve(s) %*%  gvv  %*% solve(s)  # Velicer 1974 p 565 formula (5)
		evals[,nds] <- eigen(r2) $values} 
	
} 

# mean & percentile eigenvalues for each position
means <- apply(evals, 1, mean) 
# sorting the eigenvalues for each root
for (luper in 1:Nvars) { evals[luper,] <- sort(evals[luper,]) }
percentiles <- as.matrix(evals[,round((percentile*Ndatasets)/100)])
realevals <- as.matrix(realevals)
NfactorsPA <- 0
for (root in 1:Nvars) { if (realevals[root,1] < percentiles[root,1]) { NfactorsPA <- root - 1; break } }
results <- cbind(1:Nvars,realevals,means,percentiles)
rownames(results) <- matrix((''),nrow(results),1)
colnames(results) <- c('Root', 'Real Data', 'Mean', 'Percentile')

if (verbose == TRUE) {
	
	message('\n\nPARALLEL ANALYSIS')
	if (datakind == 'correlations') message('\nThe entered data is a correlation matrix.') 
	if (randtype == 'generated')    message('\nRandomization method: generated data') 
	if (randtype == 'permuted')     message('\nRandomization method: permuted data') 
	message('\nType of correlations specified for the real data eigenvalues: ', ctype)
	# if (corkindRAND=='pearson')    ctypeRAND <- 'pearson' 
	# if (corkindRAND=='kendall')    ctypeRAND <- 'kendall' 
	# if (corkindRAND=='spearman')   ctypeRAND <- 'spearman' 
	# if (corkindRAND=='polychoric') ctypeRAND <- 'polychoric' 
	message('\nType of correlations specified for the random data eigenvalues: ', corkindRAND)
	if (corkind=='polychoric' & wholenums==0)  { 
		message('\nPolychoric correlations were specified but there are values in the')
		message('\nreal data matrix that are not whole numbers. Pearson correlations')
		message('\nwill be used instead.') 
	}
		
	if (randtype == 'permuted' & datakind == 'correlations') {
		message('\nThe permuted data option does not work when the entered data are a correlation matrix.')
		message('\nSwitching to generated random data instead of permuted random data.')
	}
	if (factormodel=='PCA')    message('\nExtraction Method: Principal Components') 
	if (factormodel=='PAF')    message('\nExtraction Method: Common Factor Analysis')
	if (factormodel=='image')  message('\nExtraction Method: Image Factor Extraction') 
	message('\nVariables = ', Nvars) 
	message('\nCases = ', Ncases) 
	message('\nNdatasets = ', Ndatasets) 
	message('\nPercentile = ', percentile) 
	message('\nCompare the Real Data eigenvalues below to the corresponding')
	message('random data Mean and Percentile eigenvalues\n')
	print(round(results,3), print.gap=3)
	message('\nThe number of factors according to the parallel analysis = ', NfactorsPA)
	message('\n')
}
rawparOutput <- list(eigenvalues=results, NfactorsPA=NfactorsPA)
return(invisible(rawparOutput))
}
