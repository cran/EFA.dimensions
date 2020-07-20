
# The salient loadings criterion for determing the number of factors
# Gorsuch, R. L. (1997a). Exploratory factor analysis: Its role in item analysis. 
#      Journal of Personality Assessment, 68, 532-560.
# numsals = The required number of salient loadings for a factor.
# salvalue = The loading value that is considered salient.
SALIENT <- function (data, salvalue=.4, numsals=3, corkind='pearson', verbose=FALSE) {
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
Nfactors <- cbind(NEVALSGT1(rdata)) # Number of eigenvalues > 1
paf.out <- PA_FA(rdata, Nfactors=Nfactors, verbose=FALSE)
loadings <- paf.out$structure
if (Nfactors > 1) loadings <- VARIMAX(loadings, verbose=FALSE)
rowmax <- cbind(apply(abs(loadings), 1, max))
nfSAL <- 0
for (lupec in 1:ncol(loadings)) {
	nsalients <- 0
	for (luper in 1:nrow(loadings)) {
		if (abs(loadings[luper,lupec]) >= salvalue & abs(loadings[luper,lupec]) == rowmax[luper,1]) { 
		nsalients <- nsalients + 1 
		}
	}
	if (nsalients >= numsals)  nfSAL <- nfSAL + 1 
}
if (verbose == TRUE) {
	message('\n\nThe salient loadings criterion for determining the number of components:')
	message('\nThe salient loading value = ', salvalue)
	message('\nThe required number salient loadings = ', numsals)
	message('\nThe loading matrix:\n')
	print(round(loadings,2))
	message('\nThe number of components according to the salient loadings criterion = ', nfSAL, '\n')
}
return(invisible(nfSAL))
}
