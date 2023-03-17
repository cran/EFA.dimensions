
# three methods of assessing the factorability of a correlation matrix or raw data set

FACTORABILITY <- function (data, corkind='pearson', Ncases=NULL, verbose=TRUE) {

data <- MISSING_DROP(data)

cnoms <- colnames(data) # get colnames

Nvars  <- ncol(data)

# set up cormat
cordat <- setupcormat(data, corkind=corkind, Ncases=Ncases)
cormat <- cordat$cormat
ctype  <- cordat$ctype
Ncases <- cordat$Ncases


# the determinant of the correlation matrix should be > 0.00001 for factorability
detcor <- det(cormat)

# Bartlett's identity matrix test
# cortest.bartlett(cormat)  # from the psych package, same results
# message('\nThe Bartlett test of whether a correlation matrix is significantly different')
# message('\nfrom an identity matrix (wherein all of the off-diagonal elements are zero).') 


chi2 <- (Ncases - 1 - (2 * Nvars + 5) / 6) * -log(detcor) 
df <- Nvars * (Nvars - 1) / 2
pvalue <- pchisq(chi2, df, lower.tail=F)

Rinv <- solve(cormat)
Rpart <- cov2cor(Rinv) 
cormat_sq <- cormat^2
Rpart_sq  <- Rpart^2

# overall KMO
KMOnum <- sum(cormat_sq) - sum(diag(cormat_sq))
KMOdenom <- KMOnum + (sum(Rpart_sq) - sum(diag(Rpart_sq))) 
KMO <- KMOnum / KMOdenom

diag(cormat_sq) <- 0
diag(Rpart_sq)  <- 0
KMOvars <- colSums(cormat_sq)/(colSums(cormat_sq) + colSums(Rpart_sq))
KMOvars <- matrix(KMOvars,length(KMOvars),1)
rownames(KMOvars) <- cnoms
colnames(KMOvars) <- 'Variable MSA'

if (verbose==TRUE) {

	message('\n\nThree methods of assessing the factorability of a correlation matrix or raw data set:')

	message('\nSpecified kind of correlations for this analysis: ', ctype)
	
	message('\n\nThe determinant of the correlation matrix should be > 0.00001 for factorability.')
	if (detcor >  0.00001) message('\nThe determinant is ',round(detcor,7),' which is > 0.00001, indicating factorability.')
	if (detcor <= 0.00001) message('\nThe determinant is ',round(detcor,7),' which is NOT > 0.00001, indicating NON factorability.')
	
	message('\n\nThe Bartlett test of whether a correlation matrix is significantly different')
	message('from an identity matrix (wherein all of the off-diagonal elements are zero):')
	message('\nchisq = ',chi2,'    df= ',df,'     p = ',pvalue)
	message('\nA significant difference is required for factorability.')
	
	message('\n\nThe Kaiser-Meyer-Olkin measure of sampling adequacy (MSA):\n')
	
	print(round(KMOvars,2))
	
	message('\nThe overall measure of sampling adequacy (MSA) = ',round(KMO,2))

	message('\nThe Kaiser & Rice (1974) interpretation guidelines for MSA values:')
	message('
	   KMO >= .9 is marvelous
	   KMO in the .80s is mertitorious
	   KMO in the .70s is middling
	   KMO in the .60s is medicore
	   KMO in the .50s is miserable
	   KMO < .5 is unacceptable')
	message('\nConsider excluding items with KMO values < .5 and then re-run the FACTORABILITY analyses.\n')
	
	message('\nThe overall KMO coefficient indicates the proportion of')
	message('variance in the variables that might be caused by underlying')
	message('factors. If the variables share common factors, then the')
	message('overall KMO coefficient should be close to 1.0. The overall')
	message('KMO indicates the extent to which there is at least one')
	message('latent factor underlying the variables. The overall KMO')
	message('index is considered particularly meaningful when the cases')
	message('to variables ratio is less than 1:5. The KMO coefficient for')
	message('a variable is a kind of summary index of how much a')
	message('variable overlaps with the other variables.\n')
}

factOutput <- list(chisq=chi2, df=df, pvalue=pvalue, Rimage=Rpart, KMO=KMO, KMOvars=KMOvars)

return(invisible(factOutput))

}
