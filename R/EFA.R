

PA_FA      <- function( ... ) { .Defunct("EFA", package="EFA.dimensions") }

MAXLIKE_FA <- function( ... ) { .Defunct("EFA", package="EFA.dimensions") }

IMAGE_FA   <- function( ... ) { .Defunct("EFA", package="EFA.dimensions") }



# MAXLIKE_FA <- function( ..., top=FALSE, shrink=1.0, textcolor=NULL )
# {
#   .Defunct("EFA", package="EFA.dimensions")
# }
# 
# IMAGE_FA <- function( ..., top=FALSE, shrink=1.0, textcolor=NULL )
# {
#   .Defunct("EFA", package="EFA.dimensions")
# }



EFA <- function (data, extraction = 'paf', corkind='pearson', Nfactors=NULL, Ncases=NULL, iterpaf=100, 
                 rotation='promax', ppower = 3, verbose=TRUE) {

data <- MISSING_DROP(data)

cnoms <- colnames(data) # get colnames

# set up cormat
cordat <- setupcormat(data, corkind=corkind, Ncases=Ncases)
cormat <- cordat$cormat
ctype  <- cordat$ctype
Ncases <- cordat$Ncases

Nvars <- dim(cormat)[1]

eigenvalues <- eigen(cormat)$values

varexplNOROT1 <- VarianceExplained(eigenvalues)

# Ratio of the 1st to the 2nd initial eigenvalues
evals12ratio <- eigenvalues[1] / eigenvalues[2]

if (is.null(Nfactors)) {		
	Nfactors <- EMPKC(cormat, Ncases=Ncases, verbose=FALSE)$NfactorsEMPKC
	NfactorsWasNull <- TRUE
} else {NfactorsWasNull <- FALSE}



# get the unrotated loadings for the specified extraction method

if (extraction == 'paf') {
	outp <- PA_FA(cormat=cormat, Nfactors=Nfactors, Ncases=Ncases, iterpaf=iterpaf)
	loadingsNOROT <- outp$loadingsNOROT
	communalities <- outp$communalities
}

if (extraction == 'ml') {
	outp <- MAXLIKE_FA(cormat=cormat, Nfactors=Nfactors, Ncases=Ncases)
	loadingsNOROT <- outp$loadingsNOROT
	communalities <- outp$communalities
}

if (extraction == 'image') {
	outp <- IMAGE_FA(cormat=cormat, Nfactors=Nfactors, Ncases=Ncases)
	loadingsNOROT <- outp$loadingsNOROT
	communalities <- outp$communalities
}

if (extraction == 'minres') {
	outp <- psych::fa(r=cormat, nfactors=Nfactors, n.obs=Ncases, fm='minres', rotate='none')
	loadingsNOROT <- outp$loadings[1:dim(outp$loadings)[1], 1:dim(outp$loadings)[2], drop=FALSE]
	communalities <- outp$communalities
}

if (extraction == 'uls')  {
	outp <- psych::fa(r=cormat, nfactors=Nfactors, n.obs=Ncases, fm='uls', rotate='none')
	loadingsNOROT <- outp$loadings[1:dim(outp$loadings)[1], 1:dim(outp$loadings)[2], drop=FALSE]
	communalities <- outp$communalities
}

if (extraction == 'ols')  {
	outp <- psych::fa(r=cormat, nfactors=Nfactors, n.obs=Ncases, fm='ols', rotate='none')
	loadingsNOROT <- outp$loadings[1:dim(outp$loadings)[1], 1:dim(outp$loadings)[2], drop=FALSE]
	communalities <- outp$communalities
}

if (extraction == 'wls')  {
	outp <- psych::fa(r=cormat, nfactors=Nfactors, n.obs=Ncases, fm='wls', rotate='none')
	loadingsNOROT <- outp$loadings[1:dim(outp$loadings)[1], 1:dim(outp$loadings)[2], drop=FALSE]
	communalities <- outp$communalities
}

if (extraction == 'gls')  {
	outp <- psych::fa(r=cormat, nfactors=Nfactors, n.obs=Ncases, fm='gls', rotate='none')
	loadingsNOROT <- outp$loadings[1:dim(outp$loadings)[1], 1:dim(outp$loadings)[2], drop=FALSE]
	communalities <- outp$communalities
}

if (extraction == 'alpha')  {
	outp <- psych::fa(r=cormat, nfactors=Nfactors, n.obs=Ncases, fm='alpha', rotate='none')
	loadingsNOROT <- outp$loadings[1:dim(outp$loadings)[1], 1:dim(outp$loadings)[2], drop=FALSE]
	communalities <- outp$communality  # an apparent name error in psych  communality rather than communalities
}

if (extraction == 'fullinfo' & cordat$datakind == "notcorrels")  {	
	# run only if data consists of integers
	# confirm that the check the fractional parts of each # in data is 0  (better than using is.integer)
	if (all(data %%1 == 0)) {				
		mirt_mod <- mirt::mirt(data=data, item_type=NULL, Nfactors)    # , itemtype = "graded"
		sum_mirt_mod <- mirt::summary(mirt_mod, rotate = 'varimax')
		loadingsNOROT <- sum_mirt_mod$rotF
		communalities <- sum_mirt_mod$h2
		# # mirt sometimes produces 0 values for factor loadings; changing to a v small #		
		# loadingsNOROT <- ifelse(loadingsNOROT == 0, .0000001, loadingsNOROT)
		
	} else {
		message('\n\n"data" does not consist of integers and so a full-information factor analysis cannot')
		message('be conducted. A paf extraction will be used instead.\n')
		extraction <- 'paf'
		outp <- PA_FA(cormat=cormat, Nfactors=Nfactors, Ncases=Ncases, iterpaf=iterpaf)
		loadingsNOROT <- outp$loadingsNOROT
		communalities <- outp$communalities
	}		
}

rownames(loadingsNOROT) <- cnoms
colnames(loadingsNOROT) <- c(paste('Factor ', 1:Nfactors, sep=''))

varexplNOROT2 <- VarianceExplained(eigenvalues, loadings=loadingsNOROT)

cormat_reprod <- loadingsNOROT %*% t(loadingsNOROT); diag(cormat_reprod) <- 1


# communalities
if (!is.matrix(communalities))  communalities <- as.matrix(communalities)
rownames(communalities) <- cnoms
if (ncol(communalities) == 1)  colnames(communalities) <- c('Communalities')
if (ncol(communalities) == 2)  colnames(communalities) <- c('Initial', 'Extraction')
uniquenesses <- 1 - communalities


# model statistics, based on Revelle

model <- cormat_reprod
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


fit_coefs <- FIT_COEFS(cormat, cormat_reprod, extraction=extraction, Ncases=Ncases, 
                       chisqMODEL=chisqMODEL, dfMODEL=dfMODEL, verbose=FALSE) 


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


efaOutput <- list(loadingsNOROT = loadingsNOROT,
                  loadingsROT = loadingsROT,
			            pattern = pattern,
			            structure = structure, 
			            phi = phi,
                  varexplNOROT1 = varexplNOROT1,
                  varexplNOROT2 = varexplNOROT2,
                  varexplROT = varexplROT,
			            evals12ratio = evals12ratio,
			            cormat_reprod = cormat_reprod,
			            fit_coefs = fit_coefs,
			            chisqMODEL = chisqMODEL, 
                  dfMODEL = dfMODEL, 
                  pvalue = pvalue,
                  chisqNULL = chisqNULL,
                  dfNULL = dfNULL,
                  communalities = communalities,
                  uniquenesses = uniquenesses)

  
if (verbose == TRUE) {
	
	message('\n\n\nExploratory Factor Analysis')
	
	message('\nThe specified kind of factor extraction method for this analysis: ', extraction)

	message('\nThe specified kind of correlations for this analysis: ', ctype)
	if (NfactorsWasNull == TRUE) {
		message('\nNfactors was not specified and so the EMPKC test was conducted to determine')
		message('the number of factors to extract: Nfactors = ', Nfactors,'\n')		
	} else if (NfactorsWasNull == FALSE) {
		message('\nThe specified number of factors to extraction = ', Nfactors,'\n')
	}

	message('\nCommunalities:\n')
	print(round(communalities,2))

	message('\n\nTotal Variance Explained (Initial Eigenvalues):\n')
	print(round(varexplNOROT1,2), print.gap=4)

	message('\nRatio of the 1st to the 2nd initial eigenvalues = ', round(evals12ratio,1))
	
	message('\n\nChi square = ',round(chisqMODEL,2),'   df = ',dfMODEL,'    p = ',round(pvalue,5))

	message('\n\nModel Fit Coefficients:')
	message('\nRMSR = ',  round(fit_coefs$RMSR,3))
	message('\nGFI = ',   round(fit_coefs$GFI,3))
	message('\nCAF = ',   round(fit_coefs$CAF,3))
	message('\nRMSEA = ', round(fit_coefs$RMSEA,3))
	message('\nTLI = ',   round(fit_coefs$TLI,3))
	message('\nCFI = ',   round(fit_coefs$CFI,3))
	message('\nMFI = ',   round(fit_coefs$MFI,3))
	message('\nBIC = ',   round(fit_coefs$BIC,3))
	message('\nAIC = ',   round(fit_coefs$AIC,3))
	message('\nCAIC = ',  round(fit_coefs$CAIC,3))
	message('\nSABIC = ', round(fit_coefs$SABIC,3))

	message('\n\nUnrotated Loadings:\n')
	print(round(loadingsNOROT,2), print.gap=3)

	message('\n\nTotal Variance Explained (Unrotated Sums of Squared Loadings):\n')
	print(round(efaOutput$varexplNOROT2,2), print.gap=4)		

	if (Nfactors == 1)  message('\nNo rotation because there is only one factor\n') 

	if (Nfactors > 1) {
		
		if (rotation == 'none')   message('\nRotation procedure:  No rotation')

		if (rotation != 'none') {

			message('\n\nThe specified rotation procedure: ', rotation)	
			
			if (!anyNA(loadingsROT)) {
				message('\n\nRotated Loadings:\n')	
				print(round(efaOutput$loadingsROT,2), print.gap=3) 
			}
				
			if (!anyNA(pattern)) { 
	
				message('\n\nPattern Matrix (standardized factor loadings):\n');      
				print(round(efaOutput$pattern,2), print.gap=3)
	
				message('\n\nStructure Matrix:\n');    
				print(round(efaOutput$structure,2), print.gap=3)
	
				message('\n\nFactor Correlations:\n'); 
				print(round(efaOutput$phi,2), print.gap=3)
			}
		
			message('\n\nTotal Variance Explained (Rotation Sums of Squared Loadings):\n')
			print(round(efaOutput$varexplROT,2), print.gap=4)		
		}
	}
}

return(invisible(efaOutput))
}






