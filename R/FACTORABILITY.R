
# three methods of assessing the factorability of a correlation matrix or raw data set
FACTORABILITY <- function (data, corkind='pearson', Ncases=NULL, verbose=TRUE) {
cnoms <- colnames(data) # get colnames
# determine whether data is a correlation matrix
if (nrow(data) == ncol(data)) {
	if (all(diag(data==1))) {datakind = 'correlations'}} else{ datakind = 'notcorrels'}
if (datakind == 'correlations') {
	cormat <- data 
	ctype <- 'from user'	
	if (is.null(Ncases)) {
		Ncases = 200
		message('\n"data" is a correlation matrix but Ncases was not specified, so Ncases was set = 200')		
	}
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
	if (corkind=='polychoric')  {cormat <- POLYCHORIC_R(data);           ctype <- 'Polychoric'}
}
#message('\nThree methods of assessing the factorability of a correlation matrix or raw data set:')
# the determinant of the correlation matrix should be > 0.00001 for factorability
detcor <- det(cormat)
# message('\nThe determinant of the correlation matrix should be > 0.00001 for factorability.')
# if (detcor > 0.00001) message('\nThe determinant is',round(detcor,7),'which is > 0.00001, indicating factorability.')
# if (detcor <= 0.00001) message('\nThe determinant is',round(detcor,7),'which is NOT > 0.00001, indicating NON factorability.')
# Bartlett's identity matrix test
# cortest.bartlett(cormat)  # from the psych package, same results
# message('\nThe Bartlett test of whether a correlation matrix is significantly different')
# message('\nfrom an identity matrix (wherein all of the off-diagonal elements are zero).') 
if (datakind == 'notcorrels') Ncases <- nrow(data)
Nvars  <- ncol(data)
chi2 <- (Ncases - 1 - (2 * Nvars + 5) / 6) * -log(detcor) 
df <- Nvars * (Nvars - 1) / 2
pvalue <- pchisq(chi2, df, lower.tail=F)
# message('\nchisq =',chi2,'    df=',df,'     p =',round(pvalue,8))
# message('\nA significant difference is required for factorability.')
# the Kaiser-Meyer-Olkin measure of sampling adequacy (MSA) -- Kaiser & Rice (1974) 
#message('\nThe Kaiser-Meyer-Olkin measure of sampling adequacy (MSA):')
# # using the KMO function from the psych package
# KMOobj <- KMO(cormat)
# message('\nOverall measure of sampling adequacy (MSA) =',round(KMOobj$MSA,2),'\n')
# itemMSA <- matrix(KMOobj$MSAi,length(KMOobj$MSAi),1)
# rownames(itemMSA) <- cnoms
# colnames(itemMSA) <- '  Item MSA'
# print(round(itemMSA,2))
Rinv <- solve(cormat)
Rpart <- cov2cor(Rinv) 
cormat_sq <- cormat^2
Rpart_sq  <- Rpart^2
# overall KMO
KMOnum <- sum(cormat_sq) - sum(diag(cormat_sq))
KMOdenom <- KMOnum + (sum(Rpart_sq) - sum(diag(Rpart_sq))) 
KMO <- KMOnum / KMOdenom
#message('\nOverall measure of sampling adequacy (MSA) =',round(KMO,2),'\n')
# variable KMOs
diag(cormat_sq) <- 0
diag(Rpart_sq)  <- 0
KMOvars <- colSums(cormat_sq)/(colSums(cormat_sq) + colSums(Rpart_sq))
KMOvars <- matrix(KMOvars,length(KMOvars),1)
rownames(KMOvars) <- cnoms
colnames(KMOvars) <- '  Variable MSA'
if (verbose==TRUE) {
	message('\n\nThree methods of assessing the factorability of a correlation matrix or raw data set:')
	
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
