
SCREE_PLOT <- function (data, corkind='pearson') {
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
eigvals <- eigen(cormat) $values
roots <- seq(1:length(eigvals))
plot(roots, eigvals, pch=15, xlab='Root', ylab='Eigenvalue', cex.lab=1.3, col='blue', 
     type='b', main='Scree Plot')
axis(1, at=roots)
}
