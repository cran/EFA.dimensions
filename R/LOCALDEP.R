
LOCALDEP <- function (data, corkind='pearson', verbose=TRUE) {
# eigenvalues & residuals correls after partialling out the first component
nvars  <- ncol(data)
cnoms <- colnames(data) # get colnames
# determine whether data is a correlation matrix
if (nrow(data) == ncol(data)) {
	if (all(diag(data==1))) {datakind = 'correlations'}} else{ datakind = 'notcorrels'}
if (datakind == 'correlations') { rdata <- data }
if (datakind == 'notcorrels') {
	ncases <- nrow(data)
	if (anyNA(data) == TRUE) {
		data <- na.omit(data)
		message('\nCases with missing values were found and removed from the data matrix.\n')
	}
	if (corkind=='pearson')     { rdata <- cor(data, method='pearson') }
	if (corkind=='kendall')     { rdata <- cor(data, method='kendall') }
	if (corkind=='spearman')    { rdata <- cor(data, method='spearman') }
	if (corkind=='polychoric')  { rdata <- POLYCHORIC_R(data) }
}
colnames(rdata) <- rownames(rdata) <- cnoms
eigs    <- eigen(as.matrix(rdata))
eigval  <- diag(eigs$values)
eigvect <- eigs$vectors
# eigval <- diag(svd(rdata) $d)
# eigvect <- svd(rdata) $u
# v <- svd(rdata) $v
evals <- cbind(1:nvars, (diag(eigval)))
dimnames(evals) <-list(rep('', dim(evals)[1]))
colnames(evals) <- c('root',' eigenvalue')
if (min(eigval) < 0) {
     message('\nThe correlation matrix is not positive definite, which means that')
     message('\nthere are eigenvalues less than zero. Expect errors.\n')
}
# residuals correls after partialling out the first component
loadings = eigvect %*% sqrt(eigval)
a = loadings[,1:1]
partcov = as.matrix(rdata - (a %*% t(a)))
d = diag ( (1 / sqrt(diag(partcov))))
pr = d %*% partcov %*% d
colnames(pr) <- cnoms
rownames(pr) <- cnoms
# numbers of residual correlations >= particular values
totalN <- Np1 <- Np2 <- Np3 <- Np4 <- Np5 <- Np6 <- Np7 <- Np8 <- 0
for (i in 1:(ncol(pr)-1)) {
	for (j in (i+1):ncol(pr)) {
	totalN = totalN + 1; 
	if (pr[i,j] >= .1)  Np1 = Np1 + 1
	if (pr[i,j] >= .2)  Np2 = Np2 + 1 
	if (pr[i,j] >= .3)  Np3 = Np3 + 1 
	if (pr[i,j] >= .4)  Np4 = Np4 + 1 
	if (pr[i,j] >= .5)  Np5 = Np5 + 1 
	if (pr[i,j] >= .6)  Np6 = Np6 + 1 
	if (pr[i,j] >= .7)  Np7 = Np7 + 1 
	if (pr[i,j] >= .8)  Np8 = Np8 + 1	
	}
}				
if (verbose == TRUE) { 
	
	message('\n\nEigenvalues & residuals correlations after partialling out the first component\n')
	
	if (datakind == 'correlations') { message('\n The entered data is a correlation matrix.') }
	
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
	message('\n\nEigenvalues:\n')
	print(round(evals,5))
	
	message('\nRatio of the 1st to the 2nd eigenvalue = ', round((eigval[1,1] / eigval[2,2]),1))
	
	message('\n\nOriginal correlation matrix:\n')
	print(round(rdata,2))
		
	message('\nResidual correlations after partialling out the first component:\n')
	print(round(pr,2))
	message('\n The # of residual correlations >= .1 is  ', Np1, '   percentage = ', round((Np1/totalN),2) , 
	        '\n The # of residual correlations >= .2 is  ', Np2, '   percentage = ', round((Np2/totalN),2) ,
	        '\n The # of residual correlations >= .3 is  ', Np3, '   percentage = ', round((Np3/totalN),2) ,
	        '\n The # of residual correlations >= .4 is  ', Np4, '   percentage = ', round((Np4/totalN),2) ,
	        '\n The # of residual correlations >= .5 is  ', Np5, '   percentage = ', round((Np5/totalN),2) ,
	        '\n The # of residual correlations >= .6 is  ', Np6, '   percentage = ', round((Np6/totalN),2) ,
	        '\n The # of residual correlations >= .7 is  ', Np7, '   percentage = ', round((Np7/totalN),2) ,
	        '\n The # of residual correlations >= .8 is  ', Np8, '   percentage = ', round((Np8/totalN),2), '\n\n')
}
localdepOutput <- list(correlations=rdata, residcor=pr)
return(invisible(localdepOutput))
}
