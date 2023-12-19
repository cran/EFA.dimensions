
# Principal Components Analysis

PCA <- function (data, corkind='pearson', Nfactors=NULL, Ncases=NULL, rotation='promax', ppower = 3, verbose=TRUE, rotate) {

# deprecated  
if (!missing(rotate))       rotation <- rotate
  
data <- MISSING_DROP(data)

cnoms <- colnames(data) # get colnames

# set up cormat
cordat <- setupcormat(data, corkind=corkind, Ncases=Ncases)
cormat <- cordat$cormat
ctype  <- cordat$ctype
Ncases <- cordat$Ncases

eigs         <- eigen(as.matrix(cormat))
eigenvalues  <- eigs$values
eigenvectors <- eigs$vectors

varexplNOROT1 <- VarianceExplained(eigenvalues)

# Ratio of the 1st to the 2nd initial eigenvalues
evals12ratio <- eigenvalues[1] / eigenvalues[2]


if (is.null(Nfactors)) {		
	Nfactors <- EMPKC(cormat, Ncases=Ncases, verbose=FALSE)$NfactorsEMPKC
	NfactorsWasNull <- TRUE
} else {NfactorsWasNull <- FALSE}


if (Nfactors == 1) {
	loadingsNOROT <- eigenvectors[,1:Nfactors] * sqrt(eigenvalues[1:Nfactors])
} else {loadingsNOROT <- eigenvectors[,1:Nfactors] %*% sqrt(diag(eigenvalues[1:Nfactors]))}
loadingsNOROT <- as.matrix(loadingsNOROT)
rownames(loadingsNOROT) <- cnoms
colnames(loadingsNOROT) <- c(paste('Factor ', 1:Nfactors, sep=''))


cormat_reprod <- loadingsNOROT %*% t(loadingsNOROT); diag(cormat_reprod) <- 1
      

# communalities
communalities <- as.matrix(diag(loadingsNOROT %*% t(loadingsNOROT))) 
communalities <- cbind(rep(1,length(communalities)), communalities)   # as.matrix(communalities) 
rownames(communalities) <- cnoms
colnames(communalities) <- c('Initial', 'Extraction')
uniquenesses <- 1 - communalities


fit_coefs <- FIT_COEFS(cormat, cormat_reprod, extraction='PCA', Ncases=Ncases, verbose=FALSE) 



# rotation   library(GPArotation)

varexplROT <- loadingsROT <- structure <- pattern <- phi <- NA

if (Nfactors > 1) {

	# orthogonal rotations
	
	if (rotation == 'varimax') 
		loadingsROT <- VARIMAX(loadingsNOROT, verbose=FALSE)$loadingsV
	
	if (rotation == 'quartimax') 
		loadingsROT <- GPArotation::quartimax(loadingsNOROT)$loadings

	if (rotation == 'bentlerT') 
		loadingsROT <- GPArotation::bentlerT(loadingsNOROT)$loadings

	if (rotation == 'equamax') 
		loadingsROT <- psych::equamax(loadingsNOROT)$loadings

	if (rotation == 'geominT') 
		loadingsROT <- GPArotation::geominT(loadingsNOROT)$loadings

	if (rotation == 'bifactorT') 
		loadingsROT <- GPArotation::bifactorT(loadingsNOROT)$loadings

	if (rotation == 'entropy') 
		loadingsROT <- GPArotation::entropy(loadingsNOROT)$loadings

	# GPArotation::parsimax(loadings)   not an exported object from 'namespace:GPArotation'


	# oblique rotations
	
	if (rotation == 'promax' | rotation == 'PROMAX') {
		promaxOutput <- PROMAX(loadingsNOROT, ppower=ppower, verbose=FALSE)
		pattern <- promaxOutput$pattern
		structure <- promaxOutput$structure
		phi <- promaxOutput$phi
	}

	if (rotation == 'quartimin') {
		outp <- GPArotation::quartimin(loadingsNOROT)
		pattern <- outp$loadings
		phi <- outp$Phi
		structure <- pattern %*% phi
	}

	if (rotation == 'oblimin') {
		outp <- GPArotation::oblimin(loadingsNOROT)
		pattern <- outp$loadings
		phi <- outp$Phi
		structure <- pattern %*% phi
	}

	if (rotation == 'oblimax') {
		outp <- GPArotation::oblimax(loadingsNOROT)
		pattern <- outp$loadings
		phi <- outp$Phi
		structure <- pattern %*% phi
	}

	if (rotation == 'simplimax') {
		outp <- GPArotation::simplimax(loadingsNOROT)
		pattern <- outp$loadings
		phi <- outp$Phi
		structure <- pattern %*% phi
	}

	if (rotation == 'bentlerQ') {
		outp <- GPArotation::bentlerQ(loadingsNOROT)
		pattern <- outp$loadings
		phi <- outp$Phi
		structure <- pattern %*% phi
	}

	if (rotation == 'geominQ') {
		outp <- GPArotation::geominQ(loadingsNOROT)
		pattern <- outp$loadings
		phi <- outp$Phi
		structure <- pattern %*% phi
	}

	if (rotation == 'bifactorQ') {
		outp <- GPArotation::bifactorQ(loadingsNOROT)
		pattern <- outp$loadings
		phi <- outp$Phi
		structure <- pattern %*% phi
	}


	if (!anyNA(loadingsROT)) varexplROT <- VarianceExplained(eigenvalues, loadings=loadingsROT)

	if (!anyNA(structure)) {			
		colnames(structure) <- colnames(phi) <- c(paste('Factor ', 1:Nfactors, sep=''))
		rownames(phi) <- c(paste('Factor ', 1:Nfactors, sep=''))		
		varexplROT <- VarianceExplained(eigenvalues, loadings=structure)
		# when factors are correlated, sums of squared loadings cannot be added to obtain a total variance, so remove them from the output
		varexplROT <- varexplROT[,1]
	}
}


pcaOutput <- list(loadingsNOROT = loadingsNOROT,
                  loadingsROT = loadingsROT,
			            pattern = pattern,
			            structure = structure, 
			            phi = phi,
                  varexplNOROT1 = varexplNOROT1,
                  varexplROT = varexplROT,
			            evals12ratio = evals12ratio,
			            cormat_reprod = cormat_reprod,
			            fit_coefs = fit_coefs,
                  communalities = communalities,
                  uniquenesses = uniquenesses)


if (verbose == TRUE) {
     message('\n\n\nPrincipal Components Analysis')

     message('\nSpecified kind of correlations for this analysis: ', ctype)

	if (NfactorsWasNull == TRUE) {
		message('\nNfactors was not specified and so the EMPKC test was conducted to determine')
		message('the number of factors to extraction: Nfactors = ', Nfactors,'\n')		
	} else if (NfactorsWasNull == FALSE) {
		message('\nThe specified number of factors to extraction = ', Nfactors,'\n')
	}

	message('\nCommunalities:\n')
	print(round(communalities,2))

	message('\n\nTotal Variance Explained (Initial Eigenvalues):\n')
	print(round(varexplNOROT1,2), print.gap=4)

	message('\nRatio of the 1st to the 2nd initial eigenvalues = ', round(evals12ratio,1))
	
	message('\nModel Fit Coefficients:')
	message('\nRMSR = ', round(fit_coefs$RMSR,3))
	message('\nGFI = ',  round(fit_coefs$GFI,3))
	message('\nCAF = ',  round(fit_coefs$CAF,3))

	message('\nUnrotated PCA Loadings:\n')
	print(round(pcaOutput$loadingsNOROT[,1:Nfactors],2), print.gap=3)

	if (Nfactors == 1)  message('\nNo rotation because there is only one component\n')

	if (Nfactors > 1) {
		
		if (rotation == 'none')   message('\nRotation procedure:  No rotation')

		if (rotation != 'none') {

			message('\n\nThe specified rotation procedure: ', rotation)	
			
			if (!anyNA(loadingsROT)) {
				message('\n\nRotated Loadings:\n')	
				print(round(pcaOutput$loadingsROT,2), print.gap=3) 
			}
				
			if (!anyNA(pattern)) { 
	
				message('\n\nPattern Matrix (standardized factor loadings):\n');      
				print(round(pcaOutput$pattern,2), print.gap=3)
	
				message('\n\nStructure Matrix:\n');    
				print(round(pcaOutput$structure,2), print.gap=3)
	
				message('\n\nFactor Correlations:\n'); 
				print(round(pcaOutput$phi,2), print.gap=3)
			}
		
			message('\n\nTotal Variance Explained (Rotation Sums of Squared Loadings):\n')
			print(round(pcaOutput$varexplROT,2), print.gap=4)		
		}
	}
}

return(invisible(pcaOutput))

}

