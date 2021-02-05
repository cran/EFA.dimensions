
ROOTFIT <- function (data, corkind='pearson', Ncases=NULL, factormodel='PAF', verbose = 'TRUE') {

Nvars <- ncol(data)

# the analyses will be run for 1-factor to (Nvars * .6)-factor solutions
# because problems tend to occur for very small factors
# Nroots <- floor(Nvars * .6)
Nroots <- Nvars

# set up cormat
cordat <- setupcormat(data, corkind=corkind, Ncases=Ncases)
cormat <- cordat$cormat
ctype  <- cordat$ctype
Ncases <- cordat$Ncases

 
if (factormodel=='PCA' | factormodel == 'IMAGE') {

	fits <- matrix(-9999,Nroots,5)

	fits[,1] <- 1:Nroots

	evals <- eigen(cormat)$values

	if (min(evals) < 0) {
	     message('\nThe correlation matrix is not positive definite, which means that')
	     message('there are eigenvalues less than zero. Expect problems.')
	}
	
	fits[,2] <- evals[1:Nroots]

	for (root in 1:Nroots) {
       
		if (factormodel == 'PCA') {
			pcaOutput <- PCA(cormat, corkind=corkind, Nfactors=root, Ncases=Ncases, 
			                 rotate='none', verbose=FALSE)
	
			fits[root,3:5] <- cbind(pcaOutput$fit_coefficients$RMSR, 
			                        pcaOutput$fit_coefficients$GFI, 
			                        pcaOutput$fit_coefficients$CAF)
		}
		if (factormodel == 'IMAGE') {
			imagefaOutput <- IMAGE_FA(cormat, corkind=corkind, Nfactors=root, Ncases=Ncases, 
			                          rotate='none', verbose=FALSE)
	
			fits[root,3:5] <- cbind(imagefaOutput$fit_coefficients$RMSR, 
			                        imagefaOutput$fit_coefficients$GFI, 
			                        imagefaOutput$fit_coefficients$CAF)
		}
		
	}
	colnames(fits) <- c('Root','Eigenvalue','RMSR','GFI','CAF')
	rownames(fits) <- matrix((''),nrow(fits),1)
}



if (factormodel == 'PAF' | factormodel == 'ML') {

	fits <- matrix(-9999,Nroots,13)

	fits[,1] <- 1:Nroots

	evals <- eigen(cormat)$values

	if (min(evals) < 0) {
	     message('\nThe correlation matrix is not positive definite, which means that')
	     message('there are eigenvalues less than zero. Expect problems.')
	}

	fits[,2] <- evals[1:Nroots]

	for (root in 1:Nroots) { 

		dof <- 0.5 * ((Nvars - root)^2 - Nvars - root) # The degrees of freedom for the model
		if (dof < 1) {
			fits <- fits[1:(root-1),]
			if (verbose == 'TRUE') {
				message('\nThe degrees of freedom for the model with ',
				    root,' factors is < 1 and so the procedure was stopped.')}
			break
		}
		
		if (factormodel == 'PAF') 
			Output <- PA_FA(data, corkind=corkind, Nfactors=root, rotate='none', 
		                    Ncases=Ncases, iterpaf=100, verbose=FALSE) 

		if (factormodel == 'ML') {
			essaye <- try(Output <- MAXLIKE_FA(data, corkind=corkind, Nfactors=root, rotate='none', 
		                                       Ncases=Ncases, verbose=FALSE), silent=TRUE)		                                       
			if (inherits(essaye, "try-error")) {
			#	message('\nfactanal produced an error when Nfactors = ', root)
				fits <- fits[1:(root-1),]
				break
			}
		}
				
		fits[root,3:13] <- cbind(Output$fit_coefficients$RMSR, 
		                         Output$fit_coefficients$GFI,
		                         Output$fit_coefficients$CAF,
		                         Output$fit_coefficients$RMSEA,
		                         Output$fit_coefficients$TLI,
		                         Output$fit_coefficients$CFI,
		                         Output$fit_coefficients$MFI,
		                         Output$fit_coefficients$BIC,
		                         Output$fit_coefficients$AIC,
		                         Output$fit_coefficients$CAIC,
		                         Output$fit_coefficients$SABIC)
	}
	colnames(fits) <- c('Root','Eigenvalue','RMSR','GFI','CAF','RMSEA','TLI','CFI',
	                    'MFI','BIC','AIC','CAIC','SABIC')
	rownames(fits) <- matrix((''),nrow(fits),1)	
}


if (verbose == 'TRUE' & factormodel=='PCA') {
	
	message('\n\nFit coefficients for N-factor solutions')

	message('\nType of correlations used in the analyses: ', ctype)

	if (factormodel=='PCA')   { message('\nExtraction Method: Principal Components') }
	if (factormodel=='PAF')   { message('\nExtraction Method: Common Factor Analysis')}
	if (factormodel=='ML')    { message('\nExtraction Method: Maximum Likelihood Estimation') } 

	message('\nThe number of cases = ', Ncases)
	message('\nThe number of variables = ', Nvars,'\n') 

	print(round(fits,3))
}


if (verbose == 'TRUE' & (factormodel == 'PAF' | factormodel == 'ML')) {

	message('\n\nFit coefficients for N-factor solutions')

	message('\nType of correlations used in the analyses: ', ctype)

	if (factormodel=='PCA')   { message('\nExtraction Method: Principal Components') }
	if (factormodel=='PAF')   { message('\nExtraction Method: Common Factor Analysis')}
	if (factormodel=='ML')    { message('\nExtraction Method: Maximum Likelihood Estimation') } 

	message('\nThe number of cases = ', Ncases)
	message('\nThe number of variables = ', Nvars,'\n') 

	print(round(fits,3))
}

return(invisible(fits))

}

