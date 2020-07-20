
# Standard Error Scree test for number of factors
# Zoski, K., & Jurs, S. (1996). An objective counterpart to the visual scree test
# for factor analysis: the standard error scree test. 
# Educational and Psychological Measurement, 56(3), 443-451.
SESCREE <- function (data, corkind='pearson', verbose=FALSE) {
# determine whether data is a correlation matrix
if (nrow(data) == ncol(data)) {
	if (max(diag(data)) == 1 & min(diag(data)) == 1) {datakind = 'correlations'}
} else{ datakind = 'notcorrels'}
if (datakind == 'correlations')  rdata <- data 
if (datakind == 'notcorrels') {
	ncases <- nrow(data)
	if (anyNA(data) == TRUE) {
		data <- na.omit(data)
		message('\n\nCases with missing values were found and removed from the data matrix.\n\n')
	}
	if (corkind=='pearson')     rdata <- cor(data, method='pearson') 
	if (corkind=='kendall')     rdata <- cor(data, method='kendall') 
	if (corkind=='spearman')    rdata <- cor(data, method='spearman') 
	if (corkind=='polychoric')  rdata <- POLYCHORIC_R(data) 
}
eigvals <- cbind(eigen(rdata) $values)
neigvals <- nrow(eigvals)
rootnum <- 1:neigvals
eigdata <-  cbind(rootnum, eigvals)
nfSCREE <- 0
arbiter <- 1 / neigvals
results <- matrix(1,(neigvals-2),3)
for (loop in 1:(neigvals-2)) {
	y <- cbind(eigdata[loop:neigvals,2])
	n <- nrow(y)
	x <- cbind(matrix(1,n,1), eigdata[loop:neigvals,1])
	b <- solve(t(x) %*% x) %*% t(x) %*% y
	pred <- x %*% b
	sse <- sum ((y - pred)^2)
	sderest <- sqrt(sse / (n - 2))
	results[loop,] <- cbind(eigvals[loop,1], b[2,1], sderest)
	if (sderest > arbiter) {nfSCREE <- nfSCREE + 1}
}
if (verbose == TRUE) {
	message('\n\nStandard Error Scree test for number of factors:')
	message('\n\nThe SE Estimate must be > the following value')
	message('for a component to be nontrivial: ', arbiter, '\n\n')
	results2 <- cbind(1:nrow(results),  results) 
	dimnames(results2) <-list(rep('', dim(results2)[1]))
	colnames(results2) <- c('Root', '   Eigenvalue', '     Slope', '    SE Estimate')
	print(round(results2,2))
	message('\nStandard Error Scree test Number of Factors = ', nfSCREE, '\n\n')
}
return(invisible(nfSCREE))
}
