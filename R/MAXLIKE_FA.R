



# Maximum likelihood factor analysis - using factanal from stats

MAXLIKE_FA <- function (data, corkind='pearson', Nfactors=NULL, Ncases=NULL,  
                        rotate='PROMAX', ppower=4, verbose=TRUE) {

data <- MISSING_DROP(data)

cnoms <- colnames(data) # get colnames

# set up cormat
cordat <- setupcormat(data, corkind=corkind, Ncases=Ncases)
cormat <- cordat$cormat
ctype  <- cordat$ctype
Ncases <- cordat$Ncases

smcINITIAL <- 1 - (1 / diag(solve(cormat)))  # initial communalities

eigenvalues <- eigen(cormat)$values

totvarexplNOROT1 <- VarianceExplained(eigenvalues)


if (is.null(Nfactors)) {		
	Nfactors <- EMPKC(cormat, Ncases=Ncases, verbose=FALSE)$NfactorsEMPKC
	NfactorsWasNull <- TRUE
} else {NfactorsWasNull <- FALSE}


Nvars <- dim(cormat)[2]
                       
     
     
                       
# factanal often generates errors
# the code below uses fa from psych when factanal produces an error
essaye1 <- try(factanalOutput <- 
              factanal(covmat = as.matrix(cormat), n.obs = Ncases, factors = Nfactors, 
                       rotation = 'none'), silent=TRUE)

if (inherits(essaye1, "try-error") == FALSE) {
	
	chisqMODEL <- unname(factanalOutput$STATISTIC)
	
	dfMODEL <- unname(factanalOutput$dof)
	
	pvalue <- unname(factanalOutput$PVAL)
	
	loadings <- factanalOutput$loadings[1:dim(factanalOutput$loadings)[1],
	                                    1:dim(factanalOutput$loadings)[2], drop=FALSE]
	
	uniquenesses <- factanalOutput$uniquenesses
}	

if (inherits(essaye1, "try-error")) {
		
	# using fa from psych if factanal produces an error
	essaye2 <- try(faOutput <- fa(cormat, Nfactors, rotate="none", fm="mle"), silent=TRUE) 

	if (inherits(essaye2, "try-error") == FALSE) {	

		chisqMODEL <- (Ncases - 1 - (2 * Nvars  +  5) / 6 - (2 * Nfactors) / 3) * faOutput$criteria[1] # from psych ?fa page

		dfMODEL <- faOutput$dof
		
		if (dfMODEL > 0) {pvalue <- pchisq(chisqMODEL, dfMODEL, lower.tail = FALSE)} else {pvalue <- NA}

		loadings <- faOutput$Structure[1:dim(faOutput$Structure)[1],
		                                    1:dim(faOutput$Structure)[2], drop=FALSE]

		uniquenesses <- faOutput$uniquenesses
}
	if (inherits(essaye2, "try-error") & verbose == TRUE) 	
	    message('\nerrors are produced an error when Nfactors = ', Nfactors, '\n')	
}
# factanal often generates errors, esp about starting values (for uniquenesses)
# the code below uses the PCA uniquenesses from the Nfactors solution as the starting values
# when factanal produces an error
# if (inherits(essaye1, "try-error")) {
	# loadingsPCA <- PCA(cormat,Nfactors=Nfactors,Ncases=Ncases,rotate='none',verbose=FALSE)$loadingsNOROT		
	# Uni <- diag(as.matrix((cormat - loadingsPCA %*% t(loadingsPCA))))  # PCA uniquenesses
	# essaye2 <- try(factanalOutput <- factanal(covmat = as.matrix(cormat), n.obs = Ncases, factors = Nfactors, 
    # if (inherits(essaye2, "try-error") & verbose == TRUE) 	
	    # message('\nfactanal produced an error when Nfactors = ', Nfactors, '\n')
# }



# the null model
Fnull <- sum(diag((cormat))) - log(det(cormat)) - Nvars  

chisqNULL <-  Fnull * ((Ncases - 1) - (2 * Nvars + 5) / 6 )

dfNULL <- Nvars * (Nvars - 1) / 2


#loadings <- factanalOutput$loadings[1:dim(factanalOutput$loadings)[1],1:dim(factanalOutput$loadings)[2], drop=FALSE]
rownames(loadings) <- cnoms
colnames(loadings) <- c(paste('Factor ', 1:Nfactors, sep=''))

totvarexplNOROT2 <- VarianceExplained(eigenvalues, loadings=loadings)

cormat_reproduced <- loadings %*% t(loadings); diag(cormat_reproduced) <- 1

fit_coefficients <- FIT_COEFS(cormat, cormat_reproduced, factormodel='ML', Ncases=Ncases, 
                              chisqMODEL=chisqMODEL, dfMODEL=dfMODEL, verbose=FALSE) 
       
communal <- as.matrix(diag(loadings %*% t(loadings))) 

Communalities <- cbind(smcINITIAL, communal) 
rownames(Communalities) <- cnoms
colnames(Communalities) <- c('Initial', 'Extraction')


# uniquenesses <- factanalOutput$uniquenesses



if (rotate=='none')  { 
	maxlikeOutput <- list(totvarexplNOROT1=totvarexplNOROT1,
	                      totvarexplNOROT2=totvarexplNOROT2,
	                      loadingsNOROT=loadings) 
}

if (rotate=='VARIMAX' | rotate=='PROMAX') {
		
	if (Nfactors==1) {
		maxlikeOutput <- list(totvarexplNOROT1=totvarexplNOROT1, 
		                      totvarexplNOROT2=totvarexplNOROT2, 
		                      loadingsNOROT=loadings, 
		                      loadingsROT=loadings, 
		                      structure=loadings, 
		                      pattern=loadings) 
	} 
	if (Nfactors > 1) {
		if (rotate=='VARIMAX') { 
			varimaxOutput <- VARIMAX(loadings,verbose=FALSE)

			totvarexplROT <- VarianceExplained(eigenvalues, loadings=varimaxOutput$loadingsV)

			maxlikeOutput <- list(totvarexplNOROT1=varimaxOutput$totvarexplNOROT1, 
		                    	  totvarexplNOROT2=totvarexplNOROT2, 
			                      totvarexplROT=totvarexplROT,
			                      loadingsNOROT=varimaxOutput$loadingsNOROT, 
			                      loadingsV=varimaxOutput$loadingsV)  
		} 
		if (rotate=='PROMAX')  { 
			promaxOutput  <- PROMAX(loadings,verbose=FALSE)

			totvarexplROT <- VarianceExplained(eigenvalues, loadings=promaxOutput$structure)
			# When factors are correlated, sums of squared loadings cannot be added to obtain a total variance, so remove them from the output
			totvarexplROT <-totvarexplROT[,1]

			maxlikeOutput <- list(totvarexplNOROT1=promaxOutput$totvarexplNOROT1, 
		                	      totvarexplNOROT2=totvarexplNOROT2, 
			                      totvarexplROT=totvarexplROT,
			                      loadingsNOROT=promaxOutput$loadingsNOROT, 
			                      pattern=promaxOutput$pattern, 
			                      structure=promaxOutput$structure, 
			                      phi=promaxOutput$phi) 
		}
	}
}

maxlikeOutput <- c(maxlikeOutput, 
		           list(cormat_reproduced=cormat_reproduced, 
		                fit_coefficients=fit_coefficients,
		                chisqMODEL=chisqMODEL, 
		                dfMODEL=dfMODEL, 
		                pvalue=pvalue,
		                chisqNULL = chisqNULL,
		                dfNULL = dfNULL,
		                Communalities = Communalities,
		                uniquenesses = uniquenesses))


if (verbose == TRUE) {
	message('\n\nMaximum likelihood factor analysis')
	message('\nSpecified kind of correlations for this analysis: ', ctype)
	if (NfactorsWasNull == TRUE) {
		message('\nNfactors was not specified and so the EMPKC test was conducted to determine')
		message('\nthe number of factors to extract: Nfactors =', Nfactors,'\n')		
	} else if (NfactorsWasNull == FALSE) {
		message('\nThe specified number of factors to extract = ', Nfactors,'\n')
	}

	message('\nCommunalities:\n')
	print(round(Communalities,2))

	message('\n\nTotal Variance Explained (Initial Eigenvalues):\n')
	print(round(totvarexplNOROT1,2), print.gap=4)

	message('\n\nChi square = ',round(chisqMODEL,2),'   df = ',dfMODEL,'    p = ',round(pvalue,5))

	message('\n\nModel Fit Coefficients:')
	message('\nRMSR = ', round(fit_coefficients$RMSR,3))
	message('\nGFI = ', round(fit_coefficients$GFI,3))
	message('\nCAF = ', round(fit_coefficients$CAF,3))
	message('\nRMSEA = ', round(fit_coefficients$RMSEA,3))
	message('\nTLI = ', round(fit_coefficients$TLI,3))
	message('\nCFI = ', round(fit_coefficients$CFI,3))
	message('\nMFI = ', round(fit_coefficients$MFI,3))
	message('\nBIC = ', round(fit_coefficients$BIC,3))
	message('\nAIC = ', round(fit_coefficients$AIC,3))
	message('\nCAIC = ', round(fit_coefficients$CAIC,3))
	message('\nSABIC = ', round(fit_coefficients$SABIC,3))
	
	message('\n\nUnrotated Maximum Likelihood Loadings:\n')
	print(round(loadings[,1:Nfactors],2), print.gap=3)

	message('\nTotal Variance Explained (Unrotated Sums of Squared Loadings):\n')
	print(round(maxlikeOutput$totvarexplNOROT2,2), print.gap=4)		

	if (Nfactors==1) { message('\nNo rotation because there is only one factor\n') }

	if (Nfactors > 1) {
		if (rotate=='none')    {message('\nRotation Procedure:  No Rotation')}

		if (rotate=='VARIMAX') {
			
				message('\nVarimax Rotated Loadings:\n'); 
				print(round(maxlikeOutput$loadingsV,2), print.gap=3) 

				message('\nTotal Variance Explained (Rotation Sums of Squared Loadings):\n')
				print(round(maxlikeOutput$totvarexplROT,2), print.gap=4)		
		}
		
		if (rotate=='PROMAX')  { 

				message('\nPromax Rotation Pattern Matrix:\n');      
				print(round(maxlikeOutput$pattern,2), print.gap=3)

				message('\nPromax Rotation Structure Matrix:\n');    
				print(round(maxlikeOutput$structure,2), print.gap=3)

				message('\nRotation Sums of Squared Loadings:\n')
				print(round(maxlikeOutput$totvarexplROT,2), print.gap=4)		
	
				message('\nPromax Rotation Factor Correlations:\n'); 
				print(round(maxlikeOutput$phi,2), print.gap=3)
		}
	}

	if (Nfactors==1) { message('\nNo rotation because there is only one factor\n') }
}

return(invisible(maxlikeOutput))
}








# # # Maximum likelihood factor analysis algorithm adapted from

# # Reyment, R., Joreskog, K., & Marcus, L. F. (1996). Applied Factor Analysis 
# # in the Natural Sciences. Cambridge, MA: Cambridge University Press. see p. 308


# MAXLIKE_FA <- function (data, corkind='pearson', Nfactors=NULL, Ncases=NULL,  
                        # rotate='PROMAX', ppower = 4, verbose=TRUE) {

# cnoms <- colnames(data) # get colnames

# # determine whether data is a correlation matrix
# if (nrow(data) == ncol(data)) {
	# if (all(diag(data==1))) {datakind = 'correlations'}} else{ datakind = 'notcorrels'}
# if (datakind == 'correlations')  {
	# cormat <- as.matrix(data)
	# ctype <- 'from user'
	# if (is.null(Ncases)) message('\nNcases must be provided when data is a correlation matrix.\n')
# }
 
# if (datakind == 'notcorrels') {
	# Ncases <- nrow(data)
	# if (anyNA(data) == TRUE) {
		# data <- na.omit(data)
		# message('\nCases with missing values were found and removed from the data matrix.\n')
	# }
	# if (corkind=='pearson')     {cormat <- cor(data, method='pearson');  ctype <- 'Pearson'}
	# if (corkind=='kendall')     {cormat <- cor(data, method='kendall');  ctype <- 'Kendall'}
	# if (corkind=='spearman')    {cormat <- cor(data, method='spearman'); ctype <- 'Spearman'} 
	# if (corkind=='polychoric')  {cormat <- POLYCHORIC_R(data);           ctype <- 'Polychoric'}
# }

# # if (is.null(Nfactors)) {		
	# Nfactors <- EMPKC(cormat, Ncases=Ncases, verbose=FALSE)$NfactorsEMPKC
	# NfactorsWasNull <- TRUE
# } else {NfactorsWasNull <- FALSE}


# tolerml=.001
# iterml=100


# Rho <- cormat
# k <- Nfactors
# p <- nrow(Rho)

# # preliminary singular value decomposition of Rho
# L <- diag(svd(Rho) $d)

# A <- svd(Rho) $u

# if (Nfactors == 1) { A1 <- A[,1:k] * sqrt(L[1:k,1:k])
# }else { A1 <- A[,1:k] %*% sqrt(L[1:k,1:k]) }   # Prin. Comp. loadings 

# Uni <- diag(diag(Rho-A1 %*% t(A1)))     # Uniqueness matrix

# Rh1 <- sqrt(solve(Uni)) %*% (Rho-Uni) %*% sqrt(solve(Uni)) # Matrix to iterate

# # First estimate of Maximum Likelihood Loadings
# L <- diag(svd(Rh1) $d)

# A <- svd(Rh1) $u

# A1 <- sqrt(Uni) %*% A[,1:k] %*% sqrt(L[1:k,1:k])

# check <- tolerml
# for (i in 1:iterml) {
   # Uni <- diag(diag(Rho-A1 %*% t(A1)))
   # Rh1 <- sqrt(solve(Uni)) %*% (Rho-Uni) %*% sqrt(solve(Uni))
   # L <- diag(svd(Rh1) $d)
   # A <- svd(Rh1) $u
   # A2 <- sqrt(Uni) %*% A[,1:k] %*% sqrt(L[1:k,1:k])
   # if (max(max(abs(A1-A2))) < check) {break}
   # A1 <- A2 
# }

# FacVar <- diag(t(A1) %*% A1)

# Com <- as.matrix(diag(A1 %*% t(A1)))  # communalities
# rownames(Com) <- cnoms
# colnames(Com) <- 'Communalities'

# Uniq <- matrix(1,p,1) - Com

# Resid <- Rho - A1 %*% t(A1)

# loadings <- as.matrix(A1)
# rownames(loadings) <- cnoms
# colnames(loadings) <-  c(paste('Factor ', 1:Nfactors, sep=''))


# # model statistics, based on Revelle
# cormat_reproduced <- loadings %*% t(loadings); diag(cormat_reproduced) <- 1
# cormat_reproduced <- cor.smooth(cormat_reproduced)
# cormat <- cor.smooth(cormat)  # ensure that cormat is positive semi-definite
# Nvars <- dim(cormat)[1]
# m.inv.r <- try(solve(cormat_reproduced,cormat),silent=TRUE)
# dfMODEL <- Nvars * (Nvars - 1) / 2 - Nvars * Nfactors + (Nfactors * (Nfactors - 1) / 2)
# objective <- sum(diag((m.inv.r))) - log(det(m.inv.r)) - Nvars 
# chisqMODEL <- objective * ((Ncases - 1) - (2 * Nvars + 5) / 6 - (2 * Nfactors) / 3) # from Tucker & from factanal
# if(!is.nan(chisqMODEL)) if (chisqMODEL < 0) {chisqMODEL <- 0}  
# if (dfMODEL > 0) {pvalue <- pchisq(chisqMODEL, dfMODEL, lower.tail = FALSE)} else {pvalue <- NA}


# fit_coefficients <- FIT_COEFS(cormat, cormat_reproduced, factormodel='ML', Ncases=Ncases, 
                              # chisqMODEL=chisqMODEL, dfMODEL=dfMODEL, verbose=FALSE) 
       

# totvarexpl <- as.matrix(diag(L))
# totvarexpl <- cbind(totvarexpl,totvarexpl/dim(loadings)[1])
# totvarexpl <- cbind(totvarexpl,cumsum(totvarexpl[,2]))
# colnames(totvarexpl) <- c('Eigenvalues','Proportion of Variance','Cumulative Variance')
# rownames(totvarexpl) <- c(paste('Factor ', 1:nrow(totvarexpl), sep=''))

# if (rotate=='none')  { 
	# maxlikeOutput <- list(totvarexpl=totvarexpl, loadingsNOROT=loadings, 
	                      # cormat_reproduced=cormat_reproduced, fit_coefficients=fit_coefficients,
	                      # chisqMODEL=chisqMODEL, dfMODEL=dfMODEL, pvalue=pvalue) 
# }

# if (rotate=='VARIMAX' | rotate=='PROMAX') {
		
	# if (Nfactors==1) {
		# maxlikeOutput <- list(totvarexpl=totvarexpl, loadingsNOROT=loadings, 
		                      # loadingsROT=loadings, structure=loadings, pattern=loadings, 
			                  # cormat_reproduced = cormat_reproduced, fit_coefficients=fit_coefficients,
	                          # chisqMODEL=chisqMODEL, dfMODEL=dfMODEL, pvalue=pvalue) 
	# } 
	# if (Nfactors > 1) {
		# if (rotate=='VARIMAX') { 
			# varimaxOutput <- VARIMAX(loadings,verbose=FALSE)
			# maxlikeOutput <- list(totvarexpl=varimaxOutput$totvarexpl, 
			                      # loadingsNOROT=varimaxOutput$loadingsNOROT, 
			                      # loadingsV=varimaxOutput$loadingsV, 
			                      # cormat_reproduced = cormat_reproduced,
	                              # chisqMODEL=chisqMODEL, dfMODEL=dfMODEL, pvalue=pvalue)  
		# } 
		# if (rotate=='PROMAX')  { 
			# promaxOutput <- PROMAX(loadings,verbose=FALSE)
			# maxlikeOutput <- list(totvarexpl=promaxOutput$totvarexpl, 
			                      # loadingsNOROT=promaxOutput$loadingsNOROT, 
			                      # pattern=promaxOutput$pattern, structure=promaxOutput$structure, 
			                      # phi=promaxOutput$phi, 
			                      # cormat_reproduced = cormat_reproduced,
	                              # chisqMODEL=chisqMODEL, dfMODEL=dfMODEL, pvalue=pvalue) 
		# }
	# }
# }


# if (verbose == TRUE) {
	# message('\n\nMaximum likelihood factor analysis')
	# message('\nSpecified kind of correlations for this analysis: ', ctype)
	# if (NfactorsWasNull == TRUE) {
		# message('\nNfactors was not specified and so the EMPKC test was conducted to determine')
		# message('\nthe number of factors to extract: Nfactors =', Nfactors,'\n')		
	# } else if (NfactorsWasNull == FALSE) {
		# message('\nThe specified number of factors to extract = ', Nfactors,'\n')
	# }

	# print(round(totvarexpl,2), print.gap=4)

	# message('\n\nChi square = ',round(chisqMODEL,2),'   df = ',dfMODEL,'    p = ',round(pvalue,5))

	# message('\n\nModel Fit Coefficients:')
	# message('\nRMSR = ', round(fit_coefficients$RMSR,3))
	# message('\nGFI = ', round(fit_coefficients$GFI,3))
	# message('\nRMSEA = ', round(fit_coefficients$RMSEA,3))
	# message('\nTLI = ', round(fit_coefficients$TLI,3))
	# message('\nCFI = ', round(fit_coefficients$CFI,3))
	# message('\nMFI = ', round(fit_coefficients$MFI,3))
	# message('\nBIC = ', round(fit_coefficients$BIC,3))
	# message('\nAIC = ', round(fit_coefficients$AIC,3))
	# message('\nCAIC = ', round(fit_coefficients$CAIC,3))
	
	# message('\n\nUnrotated Maximum Likelihood Loadings:\n')
	# print(round(loadings[,1:Nfactors],2), print.gap=3)

	# if (Nfactors==1) { message('\nNo rotation because there is only one factor\n') }

	# if (Nfactors > 1) {
		# if (rotate=='none')    {message('\nRotation Procedure:  No Rotation')}

		# if (rotate=='VARIMAX') {
			
				# message('\nVarimax Rotated Loadings:\n'); 
				# print(round(maxlikeOutput$loadingsV,2), print.gap=3) 

				# message('\nTotal Variance Explained (Rotation Sums of Squared Loadings):\n')
				# print(round(maxlikeOutput$totvarexpl,2), print.gap=4)		
		# }
		
		# if (rotate=='PROMAX')  { 

				# message('\nPromax Rotation Pattern Matrix:\n');      
				# print(round(maxlikeOutput$pattern,2), print.gap=3)

				# message('\nPromax Rotation Structure Matrix:\n');    
				# print(round(maxlikeOutput$structure,2), print.gap=3)

				# message('\nTotal Variance Explained (Rotation Sums of Squared Loadings):\n')
				# print(round(maxlikeOutput$totvarexpl,2), print.gap=4)		
	
				# message('\nPromax Rotation Factor Correlations:\n'); 
				# print(round(maxlikeOutput$phi,2), print.gap=3)
		# }
	# }

	# if (Nfactors==1) { message('\nNo rotation because there is only one factor\n') }
# }

# return(invisible(maxlikeOutput))
# }






