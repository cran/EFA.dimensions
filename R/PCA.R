
# Principal Components Analysis

PCA <- function (data, corkind='pearson', Nfactors=NULL, Ncases=NULL, rotate='PROMAX', ppower = 4, verbose=TRUE) {

cnoms <- colnames(data) # get colnames

# set up cormat
cordat <- setupcormat(data, corkind=corkind, Ncases=Ncases)
cormat <- cordat$cormat
ctype  <- cordat$ctype
Ncases <- cordat$Ncases

eigs         <- eigen(as.matrix(cormat))
eigenvalues  <- eigs$values
eigenvectors <- eigs$vectors


if (is.null(Nfactors)) {		
	Nfactors <- EMPKC(cormat, Ncases=Ncases, verbose=FALSE)$NfactorsEMPKC
	NfactorsWasNull <- TRUE
} else {NfactorsWasNull <- FALSE}


if (Nfactors == 1) {
	loadings <- eigenvectors[,1:Nfactors] * sqrt(eigenvalues[1:Nfactors])
} else {loadings <- eigenvectors[,1:Nfactors] %*% sqrt(diag(eigenvalues[1:Nfactors]))}
loadings <- as.matrix(loadings)
rownames(loadings) <- cnoms
colnames(loadings) <-  c(paste('Factor ', 1:Nfactors, sep=''))

cormat_reproduced <- loadings %*% t(loadings); diag(cormat_reproduced) <- 1

fit_coefficients <- FIT_COEFS(cormat, cormat_reproduced, factormodel='PCA', Ncases=Ncases, verbose=FALSE) 
       
eigenvar <- eigvalmat(eigenvalues)

communalities <- as.matrix(diag(loadings %*% t(loadings))) 
rownames(communalities) <- cnoms
colnames(communalities) <- 'Communalities'

uniquenesses <- 1 - communalities


if (rotate=='none')  { 
	pcaOutput <- list(eigenvar=eigenvar, loadingsNOROT=loadings) 
}

if (rotate=='PROMAX' | rotate=='VARIMAX') {
	
	if (Nfactors==1)  {
		pcaOutput <- list(eigenvar=eigenvar,
		                  loadingsNOROT=loadings, 
	                      loadingsROT=loadings,
	                      structure=loadings,
	                      pattern=loadings)  
	}
	if (Nfactors > 1) {
		if (rotate=='VARIMAX') { 
			varimaxOutput <- VARIMAX(loadings,verbose=FALSE)
			pcaOutput <- list(eigenvar=varimaxOutput$eigenvar,
			                  loadingsNOROT=varimaxOutput$loadings, 
			                  loadingsV=varimaxOutput$loadingsV) 
		}  
		if (rotate=='PROMAX')  { 
			promaxOutput <- PROMAX(loadings,verbose=FALSE)
			pcaOutput <- list(eigenvar=promaxOutput$eigenvar,
			                  loadingsNOROT=promaxOutput$loadings, 
			                  pattern=promaxOutput$pattern,
			                  structure=promaxOutput$structure, 
			                  phi=promaxOutput$phi) 
		}
	}
}

pcaOutput <- c(pcaOutput, 
		           list(cormat_reproduced=cormat_reproduced, 
		                fit_coefficients=fit_coefficients,
		                communalities = communalities,
		                uniquenesses = uniquenesses))


if (verbose == TRUE) {
     message('\n\nPrincipal Components Analysis')

     message('\nSpecified kind of correlations for this analysis: ', ctype)

	if (NfactorsWasNull == TRUE) {
		message('\nNfactors was not specified and so the EMPKC test was conducted to determine')
		message('the number of factors to extract: Nfactors = ', Nfactors,'\n')		
	} else if (NfactorsWasNull == FALSE) {
		message('\nThe specified number of factors to extract = ', Nfactors,'\n')
	}

	message('\nModel Fit Coefficients:')
	message('\nRMSR = ', round(fit_coefficients$RMSR,3))
	message('\nGFI = ',  round(fit_coefficients$GFI,3))
	message('\nCAF = ',  round(fit_coefficients$CAF,3))

	message('\n\nEigenvalues and factor proportions of variance:\n')
	print(round(eigenvar,2), print.gap=4)

	message('\nUnrotated PCA Loadings:\n')
	print(round(pcaOutput$loadingsNOROT[,1:Nfactors],2), print.gap=3)

		if (Nfactors==1) { message('\nNo rotation because there is only one component\n') }

		if (Nfactors > 1) {
			if (rotate=='none')    {message('\nRotation Procedure:  No Rotation\n')}

			if (rotate=='VARIMAX') {
				
				message('\nVarimax Rotated Loadings:\n'); 
				print(round(pcaOutput$loadingsV,2), print.gap=3) 

				message('\nEigenvalues and factor proportions of variance:\n')
				print(round(pcaOutput$eigenvar,2), print.gap=4)		
			}

			if (rotate=='PROMAX')  { 

				message('\nPromax Rotation Pattern Matrix:\n');      
				print(round(pcaOutput$pattern,2), print.gap=3)

				message('\nPromax Rotation Structure Matrix:\n');    
				print(round(pcaOutput$structure,2), print.gap=3)

				message('\nEigenvalues and factor proportions of variance:\n')
				print(round(pcaOutput$eigenvar,2), print.gap=4)		
	
				message('\nPromax Rotation Factor Correlations:\n'); 
				print(round(pcaOutput$phi,2), print.gap=3)
			}
		}
}

return(invisible(pcaOutput))

}

