
POLYCHORIC_R <- function (data, method='Revelle', verbose=TRUE){
if (is.integer(data) == FALSE) {
		if (all((data - trunc(data)) == 0) == FALSE)  {
		message('\nThe data matrix does not appear to consist of whole numbers and is therefore not appropriate
		       for the computation of polychoric correlations.')
		message('\nConsider stopping the program.\n')
		} 
}
if (anyNA(data) == TRUE) {
	data <- na.omit(data)
	message('\nCases with missing values were found and removed from the data matrix.\n')
}
# for (lupec in 1:ncol(data)) {
	# if (is.numeric(data[,lupec]) & is.integer(data[,lupec]) == FALSE) {
	# message('\nThe variables in the data matrix should be factors or integers. Numeric non-integer values\n')
	# message('have been detected, which are not appropriate for the computation of polychoric correlations.\n')
	# message('Consider stopping the program.\n') }
# }
	
	
# finding the max data value or # of levels (the max function does not work for factors)
nvalues <- apply(data, MARGIN = 2, function(x) max(x, na.rm=TRUE))
nvalues <- max(as.numeric(nvalues))
# use the polychoric function from the psych package (default)
if (nvalues < 9 & (is.null(method) | method=='Revelle')) {
	rpolysR <- suppressWarnings(psych::polychoric(data, smooth=TRUE))
	rpolys <- rpolysR$rho
	if (verbose == TRUE) {
		message('\n\nPolychoric correlations:\n')
		print(rpolys)
	}
}
if (max(nvalues) > 8) 
    {message('\nUsing the Fox polycor package because the maximum number of item categories is > 8\n')}
	
		
# use the hetcor function from the polycor package
if (method=='Fox' | max(nvalues) > 8) {
	data <- as.data.frame(data) # the data for hetcor must be a dataframe
	rpolysF <- polycor::hetcor(data)
	rpolys <- rpolysF$correlations
	if (verbose == TRUE) {
		message('\nTypes of correlations computed by hetcor:\n')
		rtypes <- rpolysF$type
		colnames(rtypes) <- rownames(rtypes) <- colnames(data)
		print(rtypes)
		message('\nPolychoric correlations:')
		print(rpolys)
	}
}
return(invisible(rpolys))
}
# using the polychor function instead of hetcor
# cnoms <- colnames(data) # get colnames
# rpolys <- matrix(-9999,ncol(data),ncol(data))
# for (i in 1:(ncol(data)-1)) {
# for (j in (i+1):ncol(data)) {
# rpolys[i,j] <- polychor(data[,i], data[,j],  ML=FALSE, std.err=FALSE, .9999) 
# rpolys[j,i] <- rpolys[i,j]
# }}
# diag(rpolys) <- 1
# if (min(eigen(rpolys) $values) < 0) { 	
	# message('\nOne or more negative eigenvalues exist in the matrix of')
	# message('\npolychoric correlations. The matrix was therefore smoothed')
	# message('\nby adding a ridge to the diagonal (see Werner & Wothke, 1993, p. 261).\n')
# # ridge approach = adding a constant to the diagonal so that
# # the smallest eval is > 0; Wothke 1993 p 261, and SAS Proc CALIS p 269
# constant  = .25
# increment = .25
# for (lupe in 1:1000) {
# rpolys2 = rpolys + diag(constant*diag(cbind(rpolys)))
# if ((min(eigen(rpolys2) $value)) > 0 & (min(eigen(rpolys2) $value)) < .001) {break}
# if ((min(eigen(rpolys2) $value)) <= 0) { constant = constant + increment}
# if ((min(eigen(rpolys2) $value)) >  0) { increment = increment / 2; constant = constant - increment}
# }
# rpolys <- rpolys2
# return(rpolys)
# } 
# colnames(rpolys) <- cnoms
# rownames(rpolys) <- cnoms
# return(invisible(rpolys))
# }}
