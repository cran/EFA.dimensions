
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
	vmaxres <- varimax(loadings)	 
	loadingsV <- vmaxres$loadings[]	
	rotmatV <- vmaxres$rotmat

	promaxres <- promax(loadingsV, m=ppower)

	bigA <- rotmatV %*% promaxres$rotmat

	phi  <- solve(t(bigA) %*% bigA)
	colnames(phi) <- rownames(phi) <- c(paste('Factor ', 1:ncol(loadingsV), sep=''))

	Pstructure <- promaxres$loadings %*% phi  # promax structure
	Ppattern <- promaxres$loadings[]  # promax loadings/pattern

	eigenvar <- eigvalmat(diag(t(Ppattern) %*% Ppattern))

	# reproduced correlation matrix
	cormat_reproduced <- Pstructure %*% t(Ppattern); diag(cormat_reproduced) <- 1
	

	if (verbose == TRUE) {
	
		# message('\nUnrotated Loadings:\n')
		# print(round(B,2))
	
		message('\nPromax Rotation Pattern Matrix:\n')
		print(round(Ppattern,2), print.gap=3)
	
		message('\n\nPromax Rotation Structure Matrix:\n')
		print(round(Pstructure,2), print.gap=3)

		message('\nEigenvalues and factor proportions of variance:\n')
		print(round(eigenvar,2), print.gap=4)		
	
		message('\nPromax Rotation Factor Correlations:\n')
		print(round(phi,2), print.gap=3)
	}

	promaxOutput <- list(loadingsNOROT=loadings, pattern=Ppattern, structure=Pstructure, 
	                     phi=phi, eigenvar=eigenvar, cormat_reproduced=cormat_reproduced)

}
return(invisible(promaxOutput))
}






# PROMAX <- function (loadings, ppower=4, verbose=TRUE) {

# # based on Marcus, 1993

# # the loadings are first rotated to an orthogonal varimax solution, then the
# # orthogonality of the factors is relaxed to better fit simple structure 
# # (Mulaik 1972 p 300; Rummel 1970 p 419)

# if (is.list(loadings) == 'TRUE') loadings <- loadings$loadings

# if (ncol(loadings) == 1)  {	
	# promaxOutput <- list(structure=loadings)
	# return(invisible(promaxOutput))
	# if (verbose == TRUE) message('\nWARNING: There was only one factor. Rotation was not performed.\n')
# }

# if (ncol(loadings) > 1) {

# #	loadings <- as.matrix(VARIMAX(loadings,verbose=FALSE))
	# B <- as.matrix(VARIMAX(loadings,verbose=FALSE))

	# Bstar <- B^ppower 

	# Tr <- solve(t(B) %*% B) %*% t(B) %*% Bstar

	# Tr <- Tr %*% sqrt(solve(diag(diag(t(Tr)%*%Tr))))  # Normalizes columns of Tr

	# Tpp <- solve(Tr)        # Tp' from definition

	# Tpp <- solve(sqrt(diag(diag(Tpp%*%t(Tpp))))) %*% Tpp # Normalizes rows of Tp'

	# Sr <- B %*% Tr           # Oblique reference structure

	# Phip <- Tpp %*% t(Tpp)     # Correlation between primary factors  # Table 6.VI

	# Sp <- B %*% t(Tpp)    # Primary Factor structure matrix Table 6.V

	# Pp <- B %*% solve(Tpp) # Primary Factor pattern matrix Table 6.IV

	# Pr <- B %*% solve(t(Tr))  # Reference pattern matrix

	# Phir <- t(Tr) %*% Tr    # Correlations between reference axes


	# colnames(Sp)   <-  c(paste('Factor ', 1:ncol(Sp),   sep=''))
	# colnames(Pp)   <-  c(paste('Factor ', 1:ncol(Pp),   sep=''))
	# colnames(Phip) <-  c(paste('Factor ', 1:ncol(Phip), sep=''))
	# rownames(Phip) <-  c(paste('Factor ', 1:ncol(Phip), sep=''))

	# if (verbose == TRUE) {
		# # message('\nUnrotated Loadings:\n')
		# # print(round(B,2))
		# message('\n\nPromax Rotation Structure Matrix:\n')
		# print(round(Sp,2))
		# message('\nPromax Rotation Pattern Matrix:\n')
		# print(round(Pp,2))
		# message('\nPromax Rotation Factor Correlations:\n')
		# print(round(Phip,2))
	# }

	# promaxOutput <- list(structure=Sp, pattern=Pp, correls=Phip)

	# return(invisible(promaxOutput))
# }
# }
