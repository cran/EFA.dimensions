
NEVALSGT1 <- function (data, corkind='pearson', verbose=FALSE) {
# Number of eigenvalues > 1
nvars  <- ncol(data)
# determine whether data is a correlation matrix
if (nrow(data) == ncol(data)) {
	if (max(diag(data)) == 1 & min(diag(data)) == 1) {datakind = 'correlations'}
} else{ datakind = 'notcorrels'}
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
eigvals <- cbind(eigen(rdata) $values)
nfnevalsgt1 <- 0
for (nev in 1:nrow(eigvals)) {if (eigvals[nev,] > 1) nfnevalsgt1 <- nfnevalsgt1 + 1}
if (verbose == TRUE) { 
	if (datakind == 'correlations') message('\n\nThe entered data is a correlation matrix.\n') 
	
	if (datakind == 'notcorrels') {
		message('\n\nNumber of cases in the data file = ', ncases)
		message('\nNumber of variables in the data file = ', nvars)
	
		# specification notices
		if (corkind=='pearson')    {message('\nCorrelations to be Analyzed: Pearson\n')}
		if (corkind=='kendall')    {message('\nCorrelations to be Analyzed: Kendall\n')}
		if (corkind=='spearman')   {message('\nCorrelations to be Analyzed: Spearman\n')}
		if (corkind=='polychoric') {message('\nCorrelations to be Analyzed: Polychoric\n')}
	}
	
	colnames(eigvals) <- 'Eigenvalues'
	rownames(eigvals) <- 1:length(eigvals)
	print(round(eigvals,5))
	
	message('\nThe number of eigenvalues greater than one = ', nfnevalsgt1, '\n')
}
return(invisible(nfnevalsgt1))
message('\n')
}
