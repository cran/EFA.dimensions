



MISSING_DROP <- function(data) {
	
	if (anyNA(data) == TRUE) {
		
		# total # of NAs
		totNAs <- sum(is.na(data))
		
		# number of rows (cases) with an NA
		nrowsNAs <- sum(apply(data, 1, anyNA))
		
		data <- na.omit(data)
		message('\nCases with missing values were found and removed from the data matrix.')
		
		message('\nThere were ', nrowsNAs, ' cases with missing values, and ', 
		        totNAs, ' missing values in total.\n')		
	}
	return(invisible(data))
}






# set up cormat

setupcormat <- function(data, corkind='pearson', Ncases=NULL) {

# determine whether data is a correlation matrix   # there is also a helper function in EFAtools = .is_cormat
if (nrow(data) == ncol(data)) {
	if (all(diag(data==1))) {datakind = 'correlations'}} else{ datakind = 'notcorrels'}
if (datakind == 'correlations')  {
	cormat <- as.matrix(data)
	ctype <- 'from user'
	if (is.null(Ncases)) message('\nNcases must be provided when data is a correlation matrix.\n')
}
 
if (datakind == 'notcorrels') {
	Ncases <- nrow(data)
	if (anyNA(data) == TRUE) {
		data <- na.omit(data)
		message('\nCases with missing values were found and removed from the data matrix.\n')
	}
	if (corkind=='pearson')     {cormat <- cor(data, method='pearson');  ctype <- 'Pearson'}
	if (corkind=='kendall')     {cormat <- cor(data, method='kendall');  ctype <- 'Kendall'}
	if (corkind=='spearman')    {cormat <- cor(data, method='spearman'); ctype <- 'Spearman'} 
	if (corkind=='polychoric')  {cormat <- POLYCHORIC_R(data);           ctype <- 'polychoric'}
	if (corkind=='gamma')       {cormat <- Rgamma(data, verbose=FALSE);  ctype <- 'Goodman-Kruskal gamma'}
}

# smooth cormat if it is not positive definite
eigenvalues <- eigen(cormat)$values
if (min(eigenvalues) <= 0)  cormat <- psych::cor.smooth(cormat)


setupcormatOutput <- list(cormat=cormat, ctype=ctype, Ncases=Ncases, datakind=datakind) 

return(invisible(setupcormatOutput))
}





# \cr\cr {Thompson, L. A. 2007. R (and S-PLUS) Manual to Accompany Agresti's Categorical Data 
	# Analysis (2002) 2nd edition. https://usermanual.wiki/Document/R2020and20SPLUS2020Manual20to20Accompany20Agrestis20Categorical20Data20Analysis.540471458#google_vignette}


# Goodman-Kruskal gamma - Laura Thompson
# https://usermanual.wiki/Document/R2020and20SPLUS2020Manual20to20Accompany20Agrestis20Categorical20Data20Analysis.540471458#google_vignette
# see p 25 of
# *** 2006 Thompson - S-PLUS (and R) Manual to Accompany Agresti's Categorical Data Analysis - Splusdiscrete2

Gamma.f <- function(x)
{
# x is a matrix of counts.  You can use output of crosstabs or xtabs.
n <- nrow(x)
m <- ncol(x)
res <- numeric((n-1)*(m-1))
for(i in 1:(n-1)) {
   for(j in 1:(m-1)) res[j+(m-1)*(i-1)] <- x[i,j]*sum(x[(i+1):n,(j+1):m])
}
C <- sum(res)
res <- numeric((n-1)*(m-1))
iter <- 0
for(i in 1:(n-1))
   for(j in 2:m) {
       iter <- iter+1; res[iter] <- x[i,j]*sum(x[(i+1):n,1:(j-1)])
   }
D <- sum(res)
gamma <- (C-D)/(C+D)
}



Rgamma  <- function (donnes, verbose=TRUE) {
	
	rgamma <- matrix(-9999,ncol(donnes),ncol(donnes))
	
	for (i in 1:(ncol(donnes)-1) ) {
		for (j in (i+1):ncol(donnes) ) {
			
			dat <- donnes[,c(i,j)]
			
			rgamma[i,j] <- rgamma[j,i] <- Gamma.f(table(dat[,1], dat[,2]))
		}
	}
	diag(rgamma) <- 1
	
	if (verbose == TRUE) {
		cat("\nGoodman-Kruskal gamma correlations (Thompson, 2006): \n\n" )
		print(round(rgamma,3))
	}
	
	return(invisible(rgamma))
}







VarianceExplained <- function(eigenvalues, loadings=NULL) {

	if (is.null(loadings)) {
		propvar <- eigenvalues / length(eigenvalues)
		cumvar  <- cumsum(propvar)
		totvarexpl <- cbind(eigenvalues, propvar, cumvar) 
	 	colnames(totvarexpl) <- c('Eigenvalues','Proportion of Variance','Cumulative Prop. Variance')
	} else {
		sumsqloads <- colSums(loadings**2)
		propvar <- sumsqloads / length(eigenvalues)
		cumvar  <- cumsum(propvar)
		totvarexpl <- cbind(sumsqloads, propvar, cumvar) 
	 	colnames(totvarexpl) <- c('Sums of Squared Loadings','Proportion of Variance','Cumulative Prop. Variance')
	}

	rownames(totvarexpl) <- c(paste('Factor ', 1:nrow(totvarexpl), sep=''))
	
	return(invisible(totvarexpl))
}





CAF_boc <- function(cormat, cormat_reproduced=NULL) {

	# from 2011 Lorenzo-Seva - The Hull method for selecting the number of common factors

	if (is.null(cormat_reproduced)) {bigR <- cormat
	} else {bigR <- cormat - cormat_reproduced;  diag(bigR) <- 1}

	# smooth bigR if it is not positive definite
	if (any(eigen(bigR, symmetric = TRUE, only.values = TRUE)$values <= 0)) 
		bigR <- suppressWarnings(psych::cor.smooth(bigR))
	
	# overall KMO
	Rinv <- solve(bigR)	
	Rpart <- cov2cor(Rinv)
	cormat_sq <- bigR^2
	Rpart_sq  <- Rpart^2
	KMOnum <- sum(cormat_sq) - sum(diag(cormat_sq))
	KMOdenom <- KMOnum + (sum(Rpart_sq) - sum(diag(Rpart_sq))) 
	KMO <- KMOnum / KMOdenom
	
	CAF <- 1 - KMO
	
	return(CAF)	
}







FIT_COEFS <- function(cormat, cormat_reproduced, extraction='PCA', Ncases=NULL, 
                      chisqMODEL=NULL, dfMODEL=NULL, verbose=TRUE) {
       
     # RMSR
     residuals <- cormat - cormat_reproduced 
     residuals.upper <- as.matrix(residuals[upper.tri(residuals, diag = FALSE)])
     mnsqdresid <- mean(residuals.upper^2) # mean of the off-diagonal squared residuals (as in Waller's MicroFact)
     RMSR <- sqrt(mean(residuals.upper^2)) # rmr is perhaps the more common term for this stat
     # no srmsr computation because it requires the SDs for the variables in the matrix

     
     # GFI (McDonald, 1999), & was also from Waller's MicroFact: 1 - mean-squared residual / mean-squared correlation
	 mnsqdcorrel <- mean(cormat[upper.tri(cormat, diag = FALSE)]^2)
     GFI <- 1 - (mnsqdresid / mnsqdcorrel)


	# CAF from Lorenzo-Seva, Timmerman, & Kiers (2011)
	 CAF <- CAF_boc(cormat, cormat_reproduced=cormat_reproduced)


	# # fit Revelle -- not using because it includes the diagonal values
    # r2 <-sum(cormat * cormat)
    # rstar <- cormat - cormat_reproduced
    # rstar2 <- sum(rstar * rstar)
    # fitRevelle <- 1- rstar2 / r2
	# factor.fit {psych}
	# The basic factor or principal components model is that a correlation or covariance matrix 
	# may be reproduced by the product of a factor loading matrix times its transpose: F'F or P'P. 
	# One simple index of fit is the 1 - sum squared residuals/sum squared original correlations. 
	# This fit index is used by VSS, ICLUST, etc.		
	# There are probably as many fit indices as there are psychometricians. This fit is a plausible 
	# estimate of the amount of reduction in a correlation matrix given a factor model. Note that 
	# it is sensitive to the size of the original correlations. That is, if the residuals are small 
	# but the original correlations are small, that is a bad fit.

	# fit.off Revelle - "how well are the off diagonal elements reproduced?"   sounds identical to Waller's GFI
	# PROBS with the commands below need fixing
	# np.obs = 10
	# #rstar.off <- sum(residual^2 * np.obs)  #weight the residuals by their sample size
	# rstar.off <- sum(rstar^2 * np.obs)  #weight the residuals by their sample size
	# r2.off <-(cormat * cormat * np.obs)   #weight the original by sample size
	# r2.off <- sum(r2.off) - tr(r2.off) 
	# fit.off <- 1-rstar.off/r2.off;  fit.off
  
    fitcoefsOutput <- list(RMSR=RMSR, GFI=GFI, CAF=CAF)


if (extraction != 'PCA' & extraction != 'pca' & extraction != 'IMAGE' | extraction != 'image') {

		Nvars <- dim(cormat)[1]
  
   		# the null model
	   	Fnull <- sum(diag((cormat))) - log(det(cormat)) - Nvars  
		chisqNULL <-  Fnull * ((Ncases - 1) - (2 * Nvars + 5) / 6 )
		dfNULL <- Nvars * (Nvars - 1) / 2

		RMSEA <- sqrt(max(((chisqMODEL - dfMODEL) / (Ncases - 1)),0) / dfMODEL)

		# TLI - Tucker-Lewis index (Tucker & Lewis, 1973) = 
        # NNFI - nonnormed fit index (Bentler & Bonett, 1980)
        t1 <- chisqNULL / dfNULL - chisqMODEL / dfMODEL
        t2 <- chisqNULL / dfNULL - 1 
        TLI <- 1
        if(t1 < 0 && t2 < 0) {TLI <- 1} else {TLI <- t1/t2} # lavaan   else {TLI <- 1}  
		NNFI <- TLI
		
		CFI <- ((chisqNULL - dfNULL) - (chisqMODEL - dfMODEL)) / (chisqNULL - dfNULL)

		# MacDonald & Marsh (1990) MFI = an absolute fit index that does not depend  
		# on comparison with another model  (T&F, 2001, p 700)
		MFI <- exp (-.5 * ( (chisqMODEL - dfMODEL) / Ncases))

		BIC <- chisqMODEL - dfMODEL * log(Ncases)

		# AIC Akaike Information Criteria (T&F, 2001, p 700)
		# not on a 0-1 scale; & the value/formula varies across software
		AIC <- chisqMODEL - 2 * dfMODEL

		# CAIC Consistent Akaike Information Criteria (T&F, 2001, p 700)
		# not on a 0-1 scale; & the value/formula varies across software
		CAIC <- chisqMODEL - (log(Ncases) + 1) * dfMODEL

		# SABIC -- Sample-Size Adjusted BIC (degree of parsimony fit index)		
		# Kenny (2020): "Like the BIC, the sample-size adjusted BIC or SABIC places a penalty 
		# for adding parameters based on sample size, but not as high a penalty as the BIC.  
		# Several recent simulation studies (Enders & Tofighi, 2008; Tofighi, & Enders, 2007) 
		SABIC <- chisqMODEL + log((Ncases+2) / 24) * (Nvars * (Nvars+1) / 2 - dfMODEL)

		# mirt:   SABIC <- (-2) * logLik + tmp*log((N+2)/24)

		fitcoefsOutput <- c(fitcoefsOutput, 
		                    list(RMSEA=RMSEA, TLI=TLI, CFI=CFI, MFI=MFI, BIC=BIC, AIC=AIC, 
		                         CAIC=CAIC, SABIC=SABIC))
}


     if(verbose == TRUE & (extraction == 'PCA' | extraction == 'pca' | 
        extraction == 'IMAGE' | extraction == 'image')) {
		
     	message('\n\nFit Coefficients:')
	
		message('\nRMSR = ', round(RMSR,2))
	
		message('\nGFI (McDonald) = ', round(GFI,2))

		message('\nCAF = ', round(CAF,2))
	}

     if(verbose == TRUE & (extraction != 'PCA' & extraction != 'pca' & 
        extraction != 'IMAGE' | extraction != 'image')) {
		
     	message('\n\nFit Coefficients:')
	
		message('\nRMSR = ', round(RMSR,2))
	
		message('\nGFI (McDonald) = ', round(GFI,2))

		message('\nCAF = ', round(CAF,2))
	
		message('\nRMSEA = ', round(RMSEA,3))
	
		message('\nTLI = ', round(TLI,2))
	
		message('\nCFI = ', round(CFI,2))
	
		message('\nMFI = ', round(MFI,2))
	
		message('\nAIC = ', round(AIC,2))
	
		message('\nCAIC = ', round(CAIC,2))

		message('\nBIC = ', round(BIC,2))
	
		message('\nSABIC = ', round(SABIC,2))
	}

return(invisible(fitcoefsOutput))    
}     







VARIMAX <- function (loadings, normalize = TRUE, verbose=TRUE) {

# uses the R built-in varimax function & provides additional output

if (is.list(loadings) == 'TRUE')  loadings <- loadings$loadings

if (ncol(loadings) == 1 & verbose==TRUE) {
	message('\nWARNING: There was only one factor. Rotation was not performed.\n')
}

if (ncol(loadings) > 1) {
	
	vmaxres <- varimax(loadings, normalize=normalize)  # from built-in stats
	 
	loadingsV <- vmaxres$loadings[]
	colnames(loadingsV) <-  c(paste('Factor ', 1:ncol(loadingsV), sep=''))

	rotmatV <- vmaxres$rotmat
	colnames(rotmatV) <- rownames(rotmatV) <- c(paste('Factor ', 1:ncol(loadingsV), sep=''))
		
	# reproduced correlation matrix
	cormat_reproduced <- loadingsV %*% t(loadingsV); diag(cormat_reproduced) <- 1
	

	if (verbose == TRUE) {

		message('\n\nVarimax Rotated Loadings:\n')
		print(round(loadingsV,2), print.gap=3)
		
		message('\n\nThe rotation matrix:\n')
		print(round(rotmatV,2), print.gap=)
	}
}

varimaxOutput <-  list(loadingsNOROT=loadings, loadingsV=loadingsV, rotmatV=rotmatV, 
                       cormat_reproduced=cormat_reproduced)  

return(invisible(varimaxOutput))

}






# Promax rotation

# from stata.com:
# The optional argument specifies the promax power. 
# Values smaller than 4 are recommended, but the choice is yours. Larger promax 
# powers simplify the loadings (generate numbers closer to zero and one) but 
# at the cost of additional correlation between factors. Choosing a value is 
# a matter of trial and error, but most sources find values in excess of 4 
# undesirable in practice. The power must be greater than 1 but is not 
# restricted to integers. 
# Promax rotation is an oblique rotation method that was developed before 
# the "analytical methods" (based on criterion optimization) became computationally 
# feasible. Promax rotation comprises an oblique Procrustean rotation of the 
# original loadings A toward the elementwise #-power of the orthogonal varimax rotation of A. 


PROMAX <- function (loadings, ppower=4, verbose=TRUE) {  
			
# uses the R built-in promax function & provides additional output

#if (is.list(loadings) == 'TRUE') loadings <- loadings$loadings

if (ncol(loadings) == 1)  {	
	promaxOutput <- list(loadingsNOROT=loadings, pattern=loadings, structure=loadings)
	return(invisible(promaxOutput))
	if (verbose == TRUE) message('\nWARNING: There was only one factor. Rotation was not performed.\n')
}

if (ncol(loadings) > 1) {

	# varimax
	vmaxres <- varimax(loadings, normalize=TRUE)   # SPSS normalizes them
	loadingsV <- vmaxres$loadings[]	
	rotmatV <- vmaxres$rotmat

	promaxres <- promax(loadingsV, m=ppower)

	bigA <- rotmatV %*% promaxres$rotmat

	phi  <- solve(t(bigA) %*% bigA)
	colnames(phi) <- rownames(phi) <- c(paste('Factor ', 1:ncol(loadingsV), sep=''))

	Pstructure <- promaxres$loadings %*% phi  # promax structure
	Ppattern   <- promaxres$loadings[]  # promax loadings/pattern

	# reproduced correlation matrix
	cormat_reproduced <- Pstructure %*% t(Ppattern); diag(cormat_reproduced) <- 1
	

	if (verbose == TRUE) {
	
		# message('\nUnrotated Loadings:\n')
		# print(round(B,2))
	
		message('\nPromax Rotation Pattern Matrix:\n')
		print(round(Ppattern,2), print.gap=3)
	
		message('\n\nPromax Rotation Structure Matrix:\n')
		print(round(Pstructure,2), print.gap=3)
	
		message('\nPromax Rotation Factor Correlations:\n')
		print(round(phi,2), print.gap=3)
	}

	promaxOutput <- list(loadingsNOROT=loadings, pattern=Ppattern, structure=Pstructure, 
	                     phi=phi, cormat_reproduced=cormat_reproduced)

}
return(invisible(promaxOutput))
}







# Image Factor Extraction (Gorsuch 1983, p 113; Velicer 1974, EPM, 34, 564)

IMAGE_FA <- function (cormat, Nfactors, Ncases) {

	smcINITIAL <- 1 - (1 / diag(solve(cormat)))  # initial communalities
	
	eigenvalues <- eigen(cormat)$values
	
	cnoms <- colnames(cormat)
		
	# factor pattern for image analysis Velicer 1974 p 565 formula (2)
	d <-  diag(1 / diag(solve(cormat)))
	gvv <- cormat + d %*% solve(cormat) %*% d - 2 * d
	s <- sqrt(d)                     #  Velicer 1974 p 565 formula (7)
	r2 <- solve(s) %*%  gvv  %*% solve(s)    #  Velicer 1974 p 565 formula (5)
	eigval <- diag(eigen(r2) $values)
	eigvect <- eigen(r2) $vectors
	l <- eigvect[,1:Nfactors]
	dd <- sqrt(eigval[1:Nfactors,1:Nfactors])
	
	loadingsNOROT <- as.matrix(s %*% l %*% dd)      #  Velicer 1974 p 565 formula (2)

	communalities <- as.matrix(diag(loadingsNOROT %*% t(loadingsNOROT))) 	
	communalities <- cbind(smcINITIAL, communalities) 
	rownames(communalities) <- cnoms
	colnames(communalities) <- c('Initial', 'Extraction')
		
	imageOutput <- list(loadingsNOROT=loadingsNOROT, communalities=communalities)
	
	return(invisible(imageOutput))
}







# Maximum likelihood factor analysis - using factanal from stats

MAXLIKE_FA <- function (cormat, Nfactors, Ncases) {

	smcINITIAL <- 1 - (1 / diag(solve(cormat)))  # initial communalities
	
	eigenvalues <- eigen(cormat)$values
	
	Nvars <- dim(cormat)[2]
	                       
	cnoms <- colnames(cormat)
		
	# factanal often generates errors
	# the code below uses fa from psych when factanal produces an error
	essaye1 <- try(factanalOutput <- 
	              factanal(covmat = as.matrix(cormat), n.obs = Ncases, factors = Nfactors, 
	                       rotation = 'none'), silent=TRUE)

	loadingsNOROT <- communalities <- NA	
	
	if (inherits(essaye1, "try-error") == FALSE) {
		
		# chisqMODEL <- unname(factanalOutput$STATISTIC)
		
		# dfMODEL <- unname(factanalOutput$dof)
		
		# pvalue <- unname(factanalOutput$PVAL)
		
		loadingsNOROT <- factanalOutput$loadings[1:dim(factanalOutput$loadings)[1],
		                                         1:dim(factanalOutput$loadings)[2], drop=FALSE]
		
		# uniquenesses <- factanalOutput$uniquenesses
	}	
	
	if (inherits(essaye1, "try-error")) {
			
		# using fa from psych if factanal produces an error
		essaye2 <- try(faOutput <- fa(cormat, Nfactors, rotate="Promax", fm="mle"), silent=TRUE) 
	
		if (inherits(essaye2, "try-error") == FALSE) {	
	
			# chisqMODEL <- (Ncases - 1 - (2 * Nvars  +  5) / 6 - (2 * Nfactors) / 3) * faOutput$criteria[1] # from psych ?fa page
	
			# dfMODEL <- faOutput$dof
			
			# if (dfMODEL > 0) {pvalue <- pchisq(chisqMODEL, dfMODEL, lower.tail = FALSE)} else {pvalue <- NA}
	
			loadingsNOROT <- faOutput$Structure[1:dim(faOutput$Structure)[1],
			                                    1:dim(faOutput$Structure)[2], drop=FALSE]
	
			# uniquenesses <- faOutput$uniquenesses
	}
		if (inherits(essaye2, "try-error")) 	
		    message('\nerrors are produced when Nfactors = ', Nfactors, '\n')	
	}
		
	# # the null model
	# Fnull <- sum(diag((cormat))) - log(det(cormat)) - Nvars  
	
	# chisqNULL <-  Fnull * ((Ncases - 1) - (2 * Nvars + 5) / 6 )
	
	# dfNULL <- Nvars * (Nvars - 1) / 2

	# if there are no errors
	if (!all(is.na(loadingsNOROT))) {				       
		communalities <- as.matrix(diag(loadingsNOROT %*% t(loadingsNOROT))) 	
		communalities <- cbind(smcINITIAL, communalities) 
		rownames(communalities) <- cnoms
		colnames(communalities) <- c('Initial', 'Extraction')
	}
	maxlikeOutput <- list(loadingsNOROT=loadingsNOROT, communalities = communalities)
	
	return(invisible(maxlikeOutput))
}








PA_FA <- function (cormat, Nfactors, Ncases, iterpaf=100) {

	# CFA / PAF  (Bernstein p 189; smc = from Bernstein p 104)
	
	Nvars <- dim(cormat)[1]
	
	eigenvalues <- eigen(cormat)$values
		
	cnoms <- colnames(cormat)
		
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
			loadingsNOROT <- eigvect[,1:Nfactors] * sqrt(eigval[1:Nfactors,1:Nfactors])
			communalities <- loadingsNOROT^2
		}else {
			loadingsNOROT <- eigvect[,1:Nfactors] %*% sqrt(eigval[1:Nfactors,1:Nfactors])
			communalities <- rowSums(loadingsNOROT^2) 
		}
		if (max(max(abs(communalities-smc))) < converge) { break }
		if (max(max(abs(communalities-smc))) >= converge  & iter < iterpaf) { smc <- communalities }
	}
	loadingsNOROT <- as.matrix(loadingsNOROT)
	
	communalities <- as.matrix(diag(loadingsNOROT %*% t(loadingsNOROT))) 	
	communalities <- cbind(smcINITIAL, communalities) 
	rownames(communalities) <- cnoms
	colnames(communalities) <- c('Initial', 'Extraction')
	
	pafOutput <- list(loadingsNOROT=loadingsNOROT, communalities=communalities) 
	
	return(invisible(pafOutput))
}












# Harris, C. W. On factors and factor scores. P 32, 363-379.

# Harris, C. W. (1962). Some Rao-Guttman relationships." Psychometrika, 27,  247-63. 

# HARRIS, CHESTER W. "Canonical Factor Models for the
# Description of Change." Problems in Measuring Change. (Edited by Chester W.
# Harris.) Madison: University of Wisconsin Press, 1963. Chapter 8, pp.
# 138-55. (a) 

# HARRIS, CHESTER W., editor. Problems in Measuring Change. Madison: University of Wisconsin Press, 1963. 259 pp. (b) 

# HARRIS, CHESTER W. "Some Recent Developments in Factor Analysis." Educational and Psychological Measurement 2 4 : 193-206; Summer 1964. 

# HARRIS, CHESTER W., and KAISER, HENRY F. "Oblique Factor Analytic Solutions by Orthogonal Transformations." Psychometrika 29: 347-62; December 1964.

# Guttman, L. (1953). Image theory for the structure of quantitative
# variates. Psychometrika 18, 277-296.









