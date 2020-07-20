
PA_FA <- function (data, corkind='pearson', Nfactors=NULL, iterpaf=100, 
                   rotate='promax', ppower=3, verbose=TRUE) {
# CFA / PAF  (Bernstein p 189; smc = from Bernstein p 104)
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
converge  <- .001
rpaf <- as.matrix(cormat)
smc <- 1 - (1 / diag(solve(rpaf)))
for (iter in 1:(iterpaf + 1)) {
	diag(rpaf) <- smc # putting smcs on the main diagonal of r
	eigval <-  diag((eigen(rpaf) $values))
	# substituting zero for negative eigenvalues
	for (luper in 1:nrow(eigval)) { if (eigval[luper,luper] < 0) { eigval[luper,luper] <- 0 }}
eigvect <- eigen(rpaf) $vectors
if (Nfactors == 1) {
	loadings <- eigvect[,1:Nfactors] * sqrt(eigval[1:Nfactors,1:Nfactors])
	communal <- loadings^2
}else {
	loadings <- eigvect[,1:Nfactors] %*% sqrt(eigval[1:Nfactors,1:Nfactors])
	communal <- rowSums(loadings^2) 
}
if (max(max(abs(communal-smc))) < converge) { break }
if (max(max(abs(communal-smc))) >= converge  & iter < iterpaf) { smc <- communal }
}
loadings <- as.matrix(loadings)
rownames(loadings) <- cnoms
colnames(loadings) <-  c(paste('  Factor ', 1:Nfactors, sep=''))
evalpaf  <-  cbind(diag(eigval))
rownames(evalpaf) <- 1:nrow(evalpaf)
colnames(evalpaf) <- 'Eigenvalues'
if (rotate=='none')  { pafOutput <- list(eigenvalues=evalpaf, loadingsNOROT=loadings) }
if (rotate=='promax' | rotate=='varimax' | rotate=='promax') {
		
	if (Nfactors==1) { 
		pafOutput <- list(eigenvalues=evalpaf, loadingsNOROT=loadings, 
		                  loadingsROT=loadings, structure=loadings, pattern=loadings) 
	} 
	if (Nfactors > 1) {
		if (rotate=='varimax') { 
			loadingsROT <- EFA.dimensions::VARIMAX(loadings,verbose=FALSE)
			pafOutput <- list(eigenvalues=evalpaf, loadingsNOROT=loadings, 
			                  loadingsROT=loadingsROT)  
		} 
		if (rotate=='promax')  { 
			loadingsROT <- EFA.dimensions::PROMAX(loadings,verbose=FALSE)
			pafOutput <- list(eigenvalues=evalpaf, structure=loadingsROT$structure, 
			                  pattern=loadingsROT$pattern, correls=loadingsROT$correls) 
		}
	}
}
if (verbose == TRUE) {
	message('\n\nPrincipal Axis Factor Analysis')
	message('\nSpecified kind of correlations for this analysis: ', ctype)
	if (NfactorsWasNull == TRUE) {
		message('\nNfactors was not specified and so the MAP test was conducted to determine')
		message('the number of factors to extract: Nfactors = ', Nfactors)		
	} else if (NfactorsWasNull == FALSE) {
		message('\nThe specified number of factors to extract = ', Nfactors)
	}
	if (max(max(abs(communal-smc))) < converge) {
		message('\nPAF converged in iterations = ', iter,'\n')	
		print(round(evalpaf,2))
		# message('\nCommunalities: \n')
		# print(round(communal,2))
		message('\nUnrotated PAF Loadings:\n')
		print(round(cbind(loadings, communal),2))
		} else { message('\nPAF did not converge in the following number of iterations:  ', (iter-1))
		}
	if (Nfactors==1) { message('\nNo rotation because there is only one factor\n') }
	if (Nfactors > 1) {
		if (rotate=='none')    {message('\nRotation Procedure:  No Rotation')}
		if (rotate=='varimax') {message('\nVarimax Rotated Loadings:\n'); print(round(loadingsROT,2))}
		if (rotate=='promax')  { 
		message('\nPromax Rotation Structure Matrix:\n');    print(round(loadingsROT$structure,2))
		message('\nPromax Rotation Pattern Matrix:\n');      print(round(loadingsROT$pattern,2))
		message('\nPromax Rotation Factor Correlations:\n'); print(round(loadingsROT$correls,2))
}}}
return(invisible(pafOutput))
}
