
# Image Factor Extraction (Gorsuch 1983, p 113; Velicer 1974, EPM, 34, 564)
IMAGE_FA <- function (data, corkind='pearson', Nfactors=NULL, rotate='promax', ppower=3, verbose=TRUE) {
cnoms <- colnames(data) # get colnames
# determine whether data is a correlation matrix
if (nrow(data) == ncol(data)) {
	if (all(diag(data==1))) {datakind = 'correlations'}} else{ datakind = 'notcorrels'}
if (datakind == 'correlations')  {
	cormat <- data 
	ctype <- 'from user'
}
 
if (datakind == 'notcorrels') {
	Ncases <- nrow(data)
	if (anyNA(data) == TRUE) {
		data <- na.omit(data)
		message('\nCases with missing values were found and removed from the data matrix.\n')
	}
	if (corkind=='pearson')     {cormat <- cor(data, method='pearson');  ctype <- 'Pearson'}
	if (corkind=='kendall')     {cormat <- cor(data, method='kendall');  ctype <- 'Kendall'}
	if (corkind=='spearman')    {cormat <- cor(data, method='spearman'); ctype <- 'Spearman'} 
	if (corkind=='polychoric')  {cormat <- POLYCHORIC_R(data);            ctype <- 'Polychoric'}
}
if (is.null(Nfactors)) {		
	nfactsMAP <- MAP(cormat, verbose=FALSE)
	Nfactors <- nfactsMAP$nfMAP
	NfactorsWasNull <- TRUE
} else {NfactorsWasNull <- FALSE}
d <-  diag(1 / diag(solve(cormat)))
gvv <- cormat + d %*% solve(cormat) %*% d - 2 * d
# factor pattern for image analysis Velicer 1974 p 565 formula (2)
s <- sqrt(d)                     #  Velicer 1974 p 565 formula (7)
r2 <- solve(s) %*%  gvv  %*% solve(s)    #  Velicer 1974 p 565 formula (5)
eigval <- diag(eigen(r2) $values)
eigvect <- eigen(r2) $vectors
l <- eigvect[,1:Nfactors]
dd <- sqrt(eigval[1:Nfactors,1:Nfactors])
evalimag  <-  cbind(diag(eigval))
rownames(evalimag) <- 1:nrow(evalimag)
colnames(evalimag) <- 'Eigenvalues'
loadings <- as.matrix(s %*% l %*% dd)      #  Velicer 1974 p 565 formula (2)
rownames(loadings) <- cnoms
colnames(loadings) <-  c(paste('  Factor ', 1:Nfactors, sep=''))
if (rotate=='none')   { imagefaOutput <- list(eigenvalues=evalimag, loadingsNOROT=loadings) }
if (rotate=='promax' | rotate=='varimax') {
	if (Nfactors==1) {
		imagefaOutput <- list(eigenvalues=evalimag, loadingsNOROT=loadings, loadingsROT=loadings, structure=loadings, pattern=loadings) 
		} 
	if (Nfactors > 1) {
		if (rotate=='varimax') { 
			loadingsROT <- EFA.dimensions::VARIMAX(loadings,verbose=FALSE)
			imagefaOutput <- list(eigenvalues=evalimag, loadingsNOROT=loadings, loadingsROT=loadingsROT) 
			} 
	if (rotate=='promax')  { 
			loadingsROT <- EFA.dimensions::PROMAX(loadings,verbose=FALSE)
			imagefaOutput <- list(eigenvalues=evalimag, structure=loadingsROT$structure, pattern=loadingsROT$pattern, correls=loadingsROT$correls) 
			}
}}
if (verbose == TRUE) {
	message('\n\nImage Factor Analysis:')
	message('\nSpecified kind of correlations for this analysis: ', ctype)
	if (NfactorsWasNull == TRUE) {
		message('\nNfactors was not specified and so the MAP test was conducted to determine')
		message('\nthe number of factors to extract: Nfactors = ', Nfactors,'\n')		
	} else if (NfactorsWasNull == FALSE) {
		message('\nThe specified number of factors to extract = ', Nfactors,'\n')
	}
	print(round(evalimag,2))
	message('\nUnrotated image Loadings:\n')
	print(round(loadings[,1:Nfactors],2))
	if (Nfactors==1) { message('\nNo rotation because there is only one factor\n') }
	if (Nfactors > 1) {
		if (rotate=='none')    {message('\nRotation Procedure:  No Rotation') }
		if (rotate=='varimax') {message('\nVarimax Rotated Loadings:\n'); print(round(loadingsROT,2)) }
		if (rotate=='promax')  { 
		message('\nPromax Rotation Structure Matrix:\n');    print(round(loadingsROT$structure,2))
		message('\nPromax Rotation Pattern Matrix:\n');      print(round(loadingsROT$pattern,2))
		message('\nPromax Rotation Factor Correlations:\n'); print(round(loadingsROT$correls,2))
}}}
return(invisible(imagefaOutput))
}
