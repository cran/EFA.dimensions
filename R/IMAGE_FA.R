

# Image Factor Extraction (Gorsuch 1983, p 113; Velicer 1974, EPM, 34, 564)

IMAGE_FA <- function (data, corkind='pearson', Nfactors=NULL, Ncases=NULL, rotate='PROMAX', ppower = 4, verbose=TRUE) {

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


# factor pattern for image analysis Velicer 1974 p 565 formula (2)
d <-  diag(1 / diag(solve(cormat)))
gvv <- cormat + d %*% solve(cormat) %*% d - 2 * d
s <- sqrt(d)                     #  Velicer 1974 p 565 formula (7)
r2 <- solve(s) %*%  gvv  %*% solve(s)    #  Velicer 1974 p 565 formula (5)
eigval <- diag(eigen(r2) $values)
eigvect <- eigen(r2) $vectors
l <- eigvect[,1:Nfactors]
dd <- sqrt(eigval[1:Nfactors,1:Nfactors])

loadings <- as.matrix(s %*% l %*% dd)      #  Velicer 1974 p 565 formula (2)
rownames(loadings) <- cnoms
colnames(loadings) <-  c(paste('Factor ', 1:Nfactors, sep=''))

totvarexplNOROT2 <- VarianceExplained(eigenvalues, loadings=loadings)

cormat_reproduced <- loadings %*% t(loadings); diag(cormat_reproduced) <- 1

fit_coefficients <- FIT_COEFS(cormat, cormat_reproduced, factormodel='IMAGE', Ncases=Ncases, verbose=FALSE) 
       
communal <- as.matrix(diag(loadings %*% t(loadings))) 

Communalities <- cbind(smcINITIAL, communal) 
rownames(Communalities) <- cnoms
colnames(Communalities) <- c('Initial', 'Extraction')

uniquenesses <- 1 - communal


if (rotate=='none')   { 
	imagefaOutput <- list(totvarexplNOROT1=totvarexplNOROT1, 
	                      totvarexplNOROT2=totvarexplNOROT2,
	                      loadingsNOROT=loadings) 
}
if (rotate=='PROMAX' | rotate=='VARIMAX') {

	if (Nfactors==1) {
		imagefaOutput <- list(totvarexplNOROT1=totvarexplNOROT1, 
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

			imagefaOutput <- list(totvarexplNOROT1=totvarexplNOROT1,
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

				imagefaOutput <- list(totvarexplNOROT1=promaxOutput$totvarexplNOROT1,
		                        	  totvarexplNOROT2=totvarexplNOROT2, 
			                          totvarexplROT=totvarexplROT,
			                          loadingsNOROT=promaxOutput$loadingsNOROT, 
				                      pattern=promaxOutput$pattern,
				                      structure=promaxOutput$structure, 
				                      phi=promaxOutput$phi) 
		}
	}
}


imagefaOutput <- c(imagefaOutput, 
		           list(cormat_reproduced=cormat_reproduced, 
		                fit_coefficients=fit_coefficients,
		                Communalities = Communalities,
		                uniquenesses = uniquenesses))


if (verbose == TRUE) {
	message('\n\nImage Factor Analysis:')
	
	message('\nSpecified kind of correlations for this analysis: ', ctype)

	if (NfactorsWasNull == TRUE) {
		message('\nNfactors was not specified and so the EMPKC test was conducted to determine')
		message('\nthe number of factors to extract: Nfactors = ', Nfactors,'\n')		
	} else if (NfactorsWasNull == FALSE) {
		message('\nThe specified number of factors to extract = ', Nfactors,'\n')
	}

	message('\nCommunalities:\n')
	print(round(Communalities,2))

	message('\n\nTotal Variance Explained (Initial Eigenvalues):\n')
	print(round(totvarexplNOROT1,2), print.gap=4)

	message('\nModel Fit Coefficients:')
	message('\nRMSR = ', round(fit_coefficients$RMSR,3))
	message('\nGFI = ', round(fit_coefficients$GFI,3))
	message('\nCAF = ', round(fit_coefficients$CAF,3))

	message('\nUnrotated image Loadings:\n')
	print(round(imagefaOutput$loadingsNOROT[,1:Nfactors],2), print.gap=3)

	message('\nTotal Variance Explained (Unrotated Sums of Squared Loadings):\n')
	print(round(imagefaOutput$totvarexplNOROT2,2), print.gap=4)		

	if (Nfactors==1) { message('\nNo rotation because there is only one factor\n') }

	if (Nfactors > 1) {

		if (rotate=='none')    {message('\nRotation Procedure:  No Rotation')}

		if (rotate=='VARIMAX') {
			
				message('\nVarimax Rotated Loadings:\n'); 
				print(round(imagefaOutput$loadingsV,2), print.gap=3) 

				message('\nTotal Variance Explained (Rotation Sums of Squared Loadings):\n')
				print(round(imagefaOutput$totvarexplROT,2), print.gap=4)		
		}
		
		if (rotate=='PROMAX')  { 

				message('\nPromax Rotation Pattern Matrix:\n');      
				print(round(imagefaOutput$pattern,2), print.gap=3)

				message('\nPromax Rotation Structure Matrix:\n');    
				print(round(imagefaOutput$structure,2), print.gap=3)

				message('\nRotation Sums of Squared Loadings:\n')
				print(round(imagefaOutput$totvarexplROT,2), print.gap=4)		
	
				message('\nPromax Rotation Factor Correlations:\n'); 
				print(round(imagefaOutput$phi,2), print.gap=3)
		}
	}
}
return(invisible(imagefaOutput))
}





# Image theory was proposed by Guttman (1953) as an alternative to the factor
# analysis model.

# A factor extraction method developed by Guttman and based on image theory.
# The common part of the variable, called the partial image, is defined as
# its linear regression on remaining variables, rather than a function of
# hypothetical factors.

# This technique in Exploratory Factor Analysis is based on the correlation
# matrix of predicted or dependent variables rather than actual variables. In
# this, we predict each variable from the others by using multiple
# regression.



# Gorsuch 1983 p 103

# Image analysis is also a principal factor variant in its usual application.
# As in the case of the principal axes with estimated communalities
# procedure, it is a principal factor variant in the sense that, after
# appropriate alterations are made to the correlation matrix, that matrix can
# be submitted to a principal factor program to find the desired factors. The
# program will then minimize the residuals of the particular matrix
# submitted. The differences in image analysis and other principal factor
# variants lie primarily in the alterations of the correlation matrix before
# the factors are extracted. Image factors can also be extracted by maximum
# likelihood procedures (Joreskog, 1969b).


# p 106

# images from the universe of variables are never factored, Only that part of
# the image of the variable is used in any given analysis that comes from the
# sample of variables included in the analysis, Since the initial theory
# considered only a population of variables and individuals to define the
# variable's image, this more restricted concept is called a partial image.
# With the number of variables usually included in a factor analysis, there
# is a distinct difference between partial image analysis and traditional
# common factor procedures such as principal axes, The traditional approach
# to common factor analysis defines the correlation between the variable and
# factor as including all of the overlap between the variable and that factor
# in addition to the overlap with the other variables, In the partial image
# approach, only that part of the variable which overlaps with other
# variables influences the correlation of the variable with the factor.
# Variance which overlaps with the factor but not with the other variables
# would be excluded from the correlation between the variable and the factor
# in the partial image approach but not in principal axes, The extent to
# which differences between principal axes and partial image analysis exist
# in any particular sample is partially a function of how the communalities
# are estimated, If overestimates of the communality are used, then the
# principal axes procedure will include specific variance in the factors. The
# partial image analysis does not allow this to occur. If the communality
# estimates for principal axes are lower than the actual communalities, then
# it will approach an image analysis solution. It is for this reason that
# Harris (1962) has called the factoring of a correlation matrix with squared
# multiple correlations in the diagonals a "backdoor image analysis," He does
# note that both the "backdoor" and true image analysis give quite similar
# results in many cases. They will, of course, be somewhat different since
# the communalities resulting from the principal axis solution are altered
# from the squared multiple correlations due to the attempt to fit the
# off-diagonal elements. This does not occur in image analysis since the
# off-diagonal elements have been altered to be consistent with the squared
# multiple correlations.



# Mulaik 2010 p 231

# Louis Guttman, who was long a critic of the indeterminacy of the common-
# factor-analysis model, provided a determinate alternative, known as "image
# analysis," which still retains many of the features of common-factor
# analysis.

# In developing image analysis Guttman relates the method to two infinitely
# large aggregates: (1) the aggregate of subjects, which is familiar to most
# statisticians as the population, and (2) the aggregate of tests (thought of
# as all the tests that logically could be constructed to measure some domain
# of characteristics of the subjects). Image analysis is concerned with both
# these aggregates. The population of subjects is used to define ultimate
# values for the correlation coefficients between tests, which are only
# estimated when one has a finite sample of subjects. The universe of tests
# is used to define ultimate values for the common and unique parts of test
# scores, which are only estimated when one works with a finite collection of
# tests.

# The distinguishing feature of image analysis lies in Guttman's explicit
# definition of the common part of a test score: "The common part of a
# variable is that part which is predictable by linear multiple correlation
# from all the other variables in the universe of variables." This common
# part of the variable is designated by Guttman as the image of the variable
# based on all the other variables in the universe of variables. The unique
# part of a variable is the part remaining which cannot be predicted by all
# the other variables in the universe of variables by multiple correlation.
# This unique part is designated as the anti-image of the variable. Note that
# the images and anti-images are defined with respect to a universe of
# variables. Obviously, in practice, we can never deal with the whole
# universe of variables corresponding to a particular domain of
# characteristics to measure. Consequently, we can never determine with
# complete accuracy, by multiple-correlation techniques, the images and
# anti-images of the variables in a finite set of variables to be analyzed.
# But the value of the image and anti-image concepts lies in their being
# ideal constructs, which, as we will see eventually, allow us to formulate a
# better model of factor analysis as well as to unify concepts of factor
# analysis with those in the area of measurement in the behavioral sciences.

# Partial-Image Analysis Although in practice we can never determine the
# total images of a set of variables by multiple correlation, we can
# nevertheless formulate completely determinate approximations to these
# images, still using multiple correlation.





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




# Guttman's Image analysis is not factor analysis in true sense ... Image
# analysis can be thought of as a specially weighted principal component
# analysis being halfway between PCA and the true common factor analysis.
# Image analysis is rarely used nowadays because the mentioned Guttman's
# demands for close-form definitions are being considered too positivistic.

# One should not confuse Guttman's Image analysis with Joreskog's Image
# Factor analysis having the fundamental theorem formula



# Harris (1962) was able to demonstrate a close relationship between Rao's
# canonical factor analysis and Guttman's (1953) "image analysis." (See
# Kaiser, 1963, for the most complete discussion of "image factor analysis.")
# In what has come to be regarded by many as the single greatest contribution
# to the study of the communality problem, Harris (1962) established the
# factoring of Guttman's image variance-covariance matrix as respectable
# psychometric practice. The relationships Harris discovered between
# canonical and image factor analyses gave a cogent answer to the question of
# when to stop extracting factors. An appropriate number of factors in both
# analyses is given by Guttman's stronger lower- bound. In addition, image
# factor analysis is a scale-free technique.
