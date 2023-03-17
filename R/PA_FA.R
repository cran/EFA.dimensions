
PA_FA <- function (data, corkind='pearson', Nfactors=NULL, Ncases=NULL, iterpaf=100, 
                   rotate='PROMAX', ppower = 4, verbose=TRUE) {

# CFA / PAF  (Bernstein p 189; smc = from Bernstein p 104)

data <- MISSING_DROP(data)

cnoms <- colnames(data) # get colnames

# set up cormat
cordat <- setupcormat(data, corkind=corkind, Ncases=Ncases)
cormat <- cordat$cormat
ctype  <- cordat$ctype
Ncases <- cordat$Ncases

Nvars <- dim(cormat)[1]

eigenvalues <- eigen(cormat)$values

totvarexplNOROT1 <- VarianceExplained(eigenvalues)

if (is.null(Nfactors)) {		
	Nfactors <- EMPKC(cormat, Ncases=Ncases, verbose=FALSE)$NfactorsEMPKC
	NfactorsWasNull <- TRUE
} else {NfactorsWasNull <- FALSE}

converge  <- .001
rpaf <- as.matrix(cormat)
smc <- 1 - (1 / diag(solve(rpaf)))
smcINITIAL <- smc  # initial communalities

for (iter in 1:(iterpaf + 1)) {
	diag(rpaf) <- smc # putting smcs on the main diagonal of r
	eigval <-  diag((eigen(rpaf) $values))
	# substituting zero for negative eigenvalues
	for (luper in 1:nrow(eigval)) { if (eigval[luper,luper] < 0) { eigval[luper,luper] <- 0 }}
	eigvect <- eigen(rpaf) $vectors
	if (Nfactors == 1) {
		loadings <- eigvect[,1:Nfactors] * sqrt(eigval[1:Nfactors,1:Nfactors])
		communal <- loadings^2
	}else {
		loadings <- eigvect[,1:Nfactors] %*% sqrt(eigval[1:Nfactors,1:Nfactors])
		communal <- rowSums(loadings^2) 
	}
	if (max(max(abs(communal-smc))) < converge) { break }
	if (max(max(abs(communal-smc))) >= converge  & iter < iterpaf) { smc <- communal }
}
loadings <- as.matrix(loadings)
rownames(loadings) <- cnoms
colnames(loadings) <-  c(paste('Factor ', 1:Nfactors, sep=''))

totvarexplNOROT2 <- VarianceExplained(eigenvalues, loadings=loadings)

cormat_reproduced <- loadings %*% t(loadings); diag(cormat_reproduced) <- 1

Communalities <- cbind(smcINITIAL, communal)   # as.matrix(communal) 
rownames(Communalities) <- cnoms
colnames(Communalities) <- c('Initial', 'Extraction')

uniquenesses <- 1 - communal


# model statistics, based on Revelle
model <- cormat_reproduced
#model <- cor.smooth(model)  #this replaces the next few lines with a slightly cleaner approach
#r <- cor.smooth(r)  #this makes sure that the correlation is positive semi-definite
m.inv.r <- try(solve(model,cormat),silent=TRUE)
dfMODEL <- Nvars * (Nvars - 1) / 2 - Nvars * Nfactors + (Nfactors * (Nfactors - 1) / 2)
objective <- sum(diag((m.inv.r))) - log(det(m.inv.r)) - Nvars 
chisqMODEL <- objective * ((Ncases - 1) - (2 * Nvars + 5) / 6 - (2 * Nfactors) / 3) # from Tucker & from factanal
if(!is.nan(chisqMODEL)) if (chisqMODEL < 0) {chisqMODEL <- 0}  
if (dfMODEL > 0) {pvalue <- pchisq(chisqMODEL, dfMODEL, lower.tail = FALSE)} else {pvalue <- NA}
    		
# the null model
Fnull <- sum(diag((cormat))) - log(det(cormat)) - Nvars  

chisqNULL <-  Fnull * ((Ncases - 1) - (2 * Nvars + 5) / 6 )

dfNULL <- Nvars * (Nvars - 1) / 2


fit_coefficients <- FIT_COEFS(cormat, cormat_reproduced, factormodel='PAF', Ncases=Ncases, 
                              chisqMODEL=chisqMODEL, dfMODEL=dfMODEL, verbose=FALSE) 


if (rotate=='none')  {
	pafOutput <- list(totvarexplNOROT1=totvarexplNOROT1, 
	                  totvarexplNOROT2=totvarexplNOROT2,
	                  loadingsNOROT=loadings) 
}
	
if (rotate=='VARIMAX' | rotate=='PROMAX') {
		
	if (Nfactors==1) {
		pafOutput <- list(totvarexplNOROT1=totvarexplNOROT1, 
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

			pafOutput <- list(totvarexplNOROT1=totvarexplNOROT1, 
		                	  totvarexplNOROT2=totvarexplNOROT2, 
			                  totvarexplROT=totvarexplROT,
			                  loadingsNOROT=varimaxOutput$loadingsNOROT, 
			                  loadingsV=varimaxOutput$loadingsV)  
		} 
		if (rotate=='PROMAX')  { 
			promaxOutput <- PROMAX(loadings,verbose=FALSE)

			totvarexplROT <- VarianceExplained(eigenvalues, loadings=promaxOutput$structure)
			# When factors are correlated, sums of squared loadings cannot be added to obtain a total variance, so remove them from the output
			totvarexplROT <-totvarexplROT[,1]

			pafOutput <- list(totvarexplNOROT1=totvarexplNOROT1, 
		                	  totvarexplNOROT2=totvarexplNOROT2, 
			                  totvarexplROT=totvarexplROT,
			                  loadingsNOROT=promaxOutput$loadingsNOROT, 
			                  pattern=promaxOutput$pattern,
			                  structure=promaxOutput$structure, 
			                  phi=promaxOutput$phi) 
		}
	}
}

pafOutput <- c(pafOutput, 
		           list(cormat_reproduced=cormat_reproduced, 
		                fit_coefficients=fit_coefficients,
		                STATISTIC=Fnull,
		                chisqMODEL=chisqMODEL, 
		                dfMODEL=dfMODEL, 
		                pvalue=pvalue,
		                chisqNULL = chisqNULL,
		                dfNULL = dfNULL,
		                Communalities = Communalities,
		                uniquenesses = uniquenesses))


if (verbose == TRUE) {
	message('\n\nPrincipal Axis Factor Analysis')
	message('\nSpecified kind of correlations for this analysis: ', ctype)
	if (NfactorsWasNull == TRUE) {
		message('\nNfactors was not specified and so the EMPKC test was conducted to determine')
		message('the number of factors to extract: Nfactors = ', Nfactors,'\n')		
	} else if (NfactorsWasNull == FALSE) {
		message('\nThe specified number of factors to extract = ', Nfactors,'\n')
	}

		message('\nCommunalities:\n')
		print(round(Communalities,2))

		message('\n\nTotal Variance Explained (Initial Eigenvalues):\n')
		print(round(totvarexplNOROT1,2), print.gap=4)

	if (max(max(abs(communal-smc))) < converge) {

		message('\nPAF converged in iterations = ', iter,'\n')	

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

		message('\n\nUnrotated PAF Loadings:\n')
		print(round(loadings,2), print.gap=3)

		} else { message('\nPAF did not converge in the following number of iterations:  ', (iter-1)) }

		message('\nTotal Variance Explained (Unrotated Sums of Squared Loadings):\n')
		print(round(pafOutput$totvarexplNOROT2,2), print.gap=4)		


	if (Nfactors==1) { message('\nNo rotation because there is only one factor\n') }

	if (Nfactors > 1) {
		if (rotate=='none')    {message('\nRotation Procedure:  No Rotation')}

		if (rotate=='VARIMAX') {
			
				message('\nVarimax Rotated Loadings:\n'); 
				print(round(pafOutput$loadingsV,2), print.gap=3) 

				message('\nTotal Variance Explained (Rotation Sums of Squared Loadings):\n')
				print(round(pafOutput$totvarexplROT,2), print.gap=4)		
		}
		
		if (rotate=='PROMAX')  { 

				message('\nPromax Rotation Pattern Matrix:\n');      
				print(round(pafOutput$pattern,2), print.gap=3)

				message('\nPromax Rotation Structure Matrix:\n');    
				print(round(pafOutput$structure,2), print.gap=3)

				message('\nRotation Sums of Squared Loadings:\n')
				print(round(pafOutput$totvarexplROT,2), print.gap=4)		
	
				message('\nPromax Rotation Factor Correlations:\n'); 
				print(round(pafOutput$phi,2), print.gap=3)
		}
	}
}

return(invisible(pafOutput))
}
