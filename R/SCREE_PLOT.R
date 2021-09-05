

SCREE_PLOT <- function (data, corkind='pearson', Ncases=NULL, verbose=TRUE) {

Nvars <- ncol(data)

# set up cormat
cordat <- setupcormat(data, corkind=corkind, Ncases=Ncases)
cormat <- cordat$cormat
ctype  <- cordat$ctype
Ncases <- cordat$Ncases


eigenvalues <- eigen(cormat)$values
roots <- seq(1:Nvars)
plot(roots, eigenvalues, pch=15, xlab='Root', ylab='Eigenvalue', cex.lab=1.3, col='blue', 
     type='b', main='Scree Plot')
axis(1, at=roots)


totvarexplNOROT <- VarianceExplained(eigenvalues)


if (verbose == TRUE) {

	message('\nScree Plot:')

	message('\nSpecified kind of correlations for this analysis: ', ctype)

	message('\n\nTotal Variance Explained (Initial Eigenvalues):\n')
	print(round(totvarexplNOROT,2), print.gap=4)

}

return(invisible(totvarexplNOROT))

}
