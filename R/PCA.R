
# Principal Components Analysis

PCA <- function (data, corkind='pearson', Nfactors=NULL, rotate='promax', ppower=3, verbose=TRUE) {

cnoms <- colnames(data) # get colnames

# determine whether data is a correlation matrix
if (nrow(data) == ncol(data)) {
	if (all(diag(data==1))) {datakind = 'correlations'}
} else{ datakind = 'notcorrels'}


if (datakind == 'correlations')  {
	cormat <- data 
	ctype <- 'from user'
}

if (datakind == 'notcorrels') {
	Ncases <- nrow(data)
	if (anyNA(data) == TRUE) {
		data <- na.omit(data)
		message('\nCases with missing values were found and removed from the data matrix.\n')
		message('\nCases with missing values were found and removed from the data matrix.\n')
	}
	if (corkind=='pearson')     {cormat <- cor(data, method='pearson');  ctype <- 'Pearson'}
	if (corkind=='kendall')     {cormat <- cor(data, method='kendall');  ctype <- 'Kendall'}
	if (corkind=='spearman')    {cormat <- cor(data, method='spearman'); ctype <- 'Spearman'} 
	if (corkind=='polychoric')  {cormat <- POLYCHORIC_R(data);           ctype <- 'Polychoric'}
}

if (is.null(Nfactors)) {		
	nfactsMAP <- MAP(cormat, verbose=FALSE)
	Nfactors <- nfactsMAP$nfMAP
	NfactorsWasNull <- TRUE
} else {NfactorsWasNull <- FALSE}

eigval <- diag(eigen(cormat) $values)
eigvect <- eigen(cormat) $vectors
if (Nfactors == 1) {loadings <- eigvect[,1:Nfactors] * sqrt(eigval[1:Nfactors,1:Nfactors])
}else {loadings <- eigvect[,1:Nfactors] %*% sqrt(eigval[1:Nfactors,1:Nfactors])}
loadings <- as.matrix(loadings)
rownames(loadings) <- cnoms
colnames(loadings) <-  c(paste('  Factor ', 1:Nfactors, sep=''))


evalpca  <-  cbind(diag(eigval))
rownames(evalpca) <- 1:nrow(evalpca)
colnames(evalpca) <- 'Eigenvalues'

if (rotate=='none')  { pcaOutput <- list(eigenvalues=evalpca, loadingsNOROT=loadings) }

if (rotate=='promax' | rotate=='varimax') {
	
	if (Nfactors==1)  pcaOutput <- list(eigenvalues=evalpca, loadingsNOROT=loadings, 
	                                    loadingsROT=loadings, structure=loadings, pattern=loadings)  

	if (Nfactors > 1) {
		if (rotate=='varimax') { 
			loadingsROT <- EFA.dimensions::VARIMAX(loadings,verbose=FALSE)
			pcaOutput <- list(eigenvalues=evalpca, loadingsNOROT=loadings, loadingsROT=loadingsROT) 
			}  
		if (rotate=='promax')  { 
			loadingsROT <- EFA.dimensions::PROMAX(loadings,verbose=FALSE)
			pcaOutput <- list(eigenvalues=evalpca, structure=loadingsROT$structure, 
			                  pattern=loadingsROT$pattern, correls=loadingsROT$correls) 
		}
	}
}

if (verbose == TRUE) {
     message('\n\nPrincipal Components Analysis')
     message('\nSpecified kind of correlations for this analysis: ', ctype)
	if (NfactorsWasNull == TRUE) {
		message('\nNfactors was not specified and so the MAP test was conducted to determine')
		message('the number of factors to extract: Nfactors = ', Nfactors,'\n')		
	} else if (NfactorsWasNull == FALSE) {
		message('\nThe specified number of factors to extract = ', Nfactors,'\n')
	}
	print(round(evalpca,2))
	message('\nUnrotated PCA Loadings:\n')
	print(round(loadings[,1:Nfactors],2))
		if (Nfactors==1) { message('\nNo rotation because there is only one component\n') }
		if (Nfactors > 1) {
			if (rotate=='none')    {message('\nRotation Procedure:  No Rotation\n')}
			if (rotate=='varimax') {message('\nVarimax Rotated Loadings:\n'); print(round(loadingsROT,2)) }
			if (rotate=='promax')  { 
				message('\nPromax Rotation Structure Matrix:\n');    print(round(loadingsROT$structure,2))
				message('\nPromax Rotation Pattern Matrix:\n');      print(round(loadingsROT$pattern,2))
				message('\nPromax Rotation Factor Correlations:\n'); print(round(loadingsROT$correls,2))
			}
		}
}

return(invisible(pcaOutput))

}

