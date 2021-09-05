

# set up cormat

setupcormat <- function(data, corkind='pearson', Ncases=NULL) {

# determine whether data is a correlation matrix
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

if (is.null(cormat_reproduced)) {bigR <- cormat
} else {bigR <- cormat - cormat_reproduced;  diag(bigR) <- 1}

	Rinv <- solve(bigR)
	Rpart <- cov2cor(Rinv) 
	cormat_sq <- bigR^2
	Rpart_sq  <- Rpart^2
	
	# overall KMO
	KMOnum <- sum(cormat_sq) - sum(diag(cormat_sq))
	KMOdenom <- KMOnum + (sum(Rpart_sq) - sum(diag(Rpart_sq))) 
	KMO <- KMOnum / KMOdenom
	
	CAF <- 1 - KMO
	
	return(CAF)
	
}






FIT_COEFS <- function(cormat, cormat_reproduced, factormodel='PCA', Ncases=NULL, 
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


if(factormodel == 'PAF' | factormodel == 'ML') {

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


     if(verbose == TRUE & (factormodel == 'PCA'| factormodel == 'IMAGE')) {
		
     	message('\n\nFit Coefficients:')
	
		message('\nRMSR = ', round(RMSR,2))
	
		message('\nGFI (McDonald) = ', round(GFI,2))

		message('\nCAF = ', round(CAF,2))
	}

     if(verbose == TRUE & (factormodel == 'PAF' | factormodel == 'ML')) {
		
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





