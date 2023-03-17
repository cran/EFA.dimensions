

MAP <- function (data, corkind='pearson', Ncases=NULL, verbose=TRUE) {
# Velicer's MAP test -- takes raw data or a correlation matrix

data <- MISSING_DROP(data)

data <- as.matrix(data)
     
Nvars  <- ncol(data)

# set up cormat
cordat <- setupcormat(data, corkind=corkind, Ncases=Ncases)
cormat <- cordat$cormat
ctype  <- cordat$ctype
Ncases <- cordat$Ncases
datakind <- cordat$datakind


eigs         <- eigen(as.matrix(cormat))
eigenvalues  <- eigs$values
eigenvectors <- eigs$vectors

totvarexplNOROT <- VarianceExplained(eigenvalues)

loadings <- eigenvectors %*% sqrt(diag(eigenvalues))

fmfm4  <- matrix(NA,Nvars,3)
fmfm4[,1]  <- 0:(Nvars-1)
fmfm4[1,2] <- (sum(sum(cormat^2))-Nvars)/(Nvars*(Nvars-1))
fmfm4[1,3] <- (sum(sum(cormat^4))-Nvars)/(Nvars*(Nvars-1))


# pb <- utils::txtProgressBar(min = 0, max = (Nvars - 1), style = 3) # create progress bar
for (m in 1:(Nvars - 1)) {
	#     Sys.sleep(0.1) # for the progress bar
	a <- loadings[,1:m]
	partcov <- as.matrix(cormat - tcrossprod(a,a))  # faster than as.matrix(cormat - (a %*% t(a)))
	
	if (max(partcov) > .0001) {
		d <- diag ( (1 / sqrt(diag(partcov))))
		pr <- d %*% (partcov %*% d)  # faster than d %*% partcov %*% d
		fmfm4[m+1,2] <- (sum(sum(pr^2))-Nvars)/(Nvars*(Nvars-1))
		fmfm4[m+1,3] <- (sum(sum(pr^4))-Nvars)/(Nvars*(Nvars-1))
	}	else {break}	
	
	#	rm(a,partcov,d,pr) # remove large matrices to free up memory re: R creates duplicates
	#     utils::setTxtProgressBar(pb, m) # update progress bar
}
# close(pb)


# identifying the smallest fm values & their locations
NfactorsMAP   <- which.min(na.omit(fmfm4[,2])) - 1
NfactorsMAP4  <- which.min(na.omit(fmfm4[,3])) - 1

dimnames(fmfm4) <-list(rep('', dim(fmfm4)[1]))
colnames(fmfm4) <- c('root','Avg.Corr.Sq.','Avg.Corr.power4')


if (verbose == TRUE) { 

	message("\n\nMINIMUM AVERAGE PARTIAL (MAP) TEST")
	
	if (datakind == 'correlations') message('\nThe entered data is a correlation matrix.') 
	
	message('\nNumber of cases = ', Ncases)
	
	message('\nNumber of variables = ', Nvars)
	
	message('\nSpecified kind of correlations for this analysis: ', ctype)
	
	message('\n\nTotal Variance Explained (Initial Eigenvalues):\n')
	print(round(totvarexplNOROT,2), print.gap=4)
	
	message("\nVelicer's Average Squared Correlations\n")
	print(round(fmfm4,5), print.gap=3)
	
	message('\nThe smallest average squared correlation is ', round(min(na.omit(fmfm4[,2])),5))
	message('\nThe smallest average 4rth power correlation is ', round(min(na.omit(fmfm4[,3])),5))
	
	message('\nThe number of components according to the original (1976) MAP Test is = ', NfactorsMAP,  labels = NULL)
	message('\nThe number of components according to the revised (2000) MAP Test is = ', NfactorsMAP4, labels = NULL, '\n')
}

mapOutput <- list(totvarexplNOROT=totvarexplNOROT, avgsqrs=fmfm4, NfactorsMAP=NfactorsMAP, NfactorsMAP4=NfactorsMAP4)

return(invisible(mapOutput))

message('\n')
}
