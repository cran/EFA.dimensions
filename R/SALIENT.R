
# The salient loadings criterion for determing the number of factors
# Gorsuch, R. L. (1997a). Exploratory factor analysis: Its role in item analysis. 
#      Journal of Personality Assessment, 68, 532-560.
# numsals = The required number of salient loadings for a factor.
# salvalue = The loading value that is considered salient.

SALIENT <- function (data, salvalue=.4, numsals=3, corkind='pearson', Ncases=NULL, verbose=TRUE) {

data <- MISSING_DROP(data)

# set up cormat
cordat <- setupcormat(data, corkind=corkind, Ncases=Ncases)
cormat <- cordat$cormat
ctype  <- cordat$ctype
Ncases <- cordat$Ncases


nevalsgt1Output <- NEVALSGT1(cormat, Ncases=Ncases)

Nfactors <-nevalsgt1Output$NfactorsNEVALSGT1 # Number of eigenvalues > 1

loadings <- PA_FA(cormat, Nfactors=Nfactors, Ncases=Ncases, rotate='none', verbose=FALSE)$loadingsNOROT

if (Nfactors > 1) loadings <- VARIMAX(loadings, verbose=FALSE)$loadingsV
rowmax <- cbind(apply(abs(loadings), 1, max))

NfactorsSALIENT <- 0
for (lupec in 1:ncol(loadings)) {
	nsalients <- 0
	for (luper in 1:nrow(loadings)) {
		if (abs(loadings[luper,lupec]) >= salvalue & abs(loadings[luper,lupec]) == rowmax[luper,1]) { 
		nsalients <- nsalients + 1 
		}
	}
	if (nsalients >= numsals)  NfactorsSALIENT <- NfactorsSALIENT + 1 
}

totvarexplNOROT <- nevalsgt1Output$totvarexplNOROT

if (verbose == TRUE) {
	message('\n\nNUMBER OF SALIENT LOADINGS:')
	message('\nSpecified kind of correlations for this analysis: ', ctype)
	message('\nThe salient loading value = ', salvalue)
	message('\nThe required number salient loadings = ', numsals)
	message('\nThe loading matrix:\n')
	print(round(loadings,2), print.gap=3)
	message('\nThe number of components according to the salient loadings criterion = ', NfactorsSALIENT, '\n')
}

salientOutput <- list(NfactorsSALIENT=NfactorsSALIENT, totvarexplNOROT=totvarexplNOROT)

return(invisible(salientOutput))

}
