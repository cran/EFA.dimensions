
LOCALDEP <- function (data, corkind='pearson', Ncases=NULL, verbose=TRUE) {
# eigenvalues & residuals correls after partialling out the first component

Nvars  <- ncol(data)

cnoms <- colnames(data) # get colnames

# set up cormat
cordat <- setupcormat(data, corkind=corkind, Ncases=Ncases)
cormat <- cordat$cormat
ctype  <- cordat$ctype
Ncases <- cordat$Ncases
datakind <- cordat$datakind


colnames(cormat) <- rownames(cormat) <- cnoms

eigs         <- eigen(as.matrix(cormat))
eigenvalues  <- eigs$values
eigenvectors <- eigs$vectors

totvarexplNOROT <- VarianceExplained(eigenvalues)

# residuals correls after partialling out the first component
loadings = eigenvectors %*% sqrt(diag(eigenvalues))
a = loadings[,1:1]
partcov = as.matrix(cormat - (a %*% t(a)))
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
	
	if (datakind == 'correlations') { message('\nThe entered data is a correlation matrix') }
	
	if (datakind == 'notcorrels') {
		message('\nNumber of cases = ', Ncases)
		message('\nNumber of variables = ', Nvars)
		# message('\nSummary statistics for the data file:\n')
		# print(summary(data))
	}

	message('\nSpecified kind of correlations for this analysis: ', ctype)

	message('\n\nTotal Variance Explained (Initial Eigenvalues):\n')
	print(round(totvarexplNOROT,2), print.gap=4)
	
	message('\nRatio of the 1st to the 2nd eigenvalue = ', round((eigenvalues[1] / eigenvalues[2]),1))
	
	message('\n\nOriginal correlation matrix:\n')
	print(round(cormat,2))
		
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
localdepOutput <- list(correlations=cormat, residcor=pr, eigenvalues=eigenvalues)
return(invisible(localdepOutput))
}
