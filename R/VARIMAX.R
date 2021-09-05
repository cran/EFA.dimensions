

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




# # VARIMAX <- function (loadings, verbose=TRUE) {
# # VARIMAX
# # Marcus, 1993, in Reyment & Joreskog, 
# # Applied Factor Analysis in the Natural Sciences, CUP.
# # This procedure follows algorithm as spelled out in Harman (1960) in Chapter 14, section 4.  
# if (is.list(loadings) == 'TRUE')  loadings <- loadings$loadings
# if (ncol(loadings) == 1 & verbose==TRUE) {
	# message('\nWARNING: There was only one factor. Rotation was not performed.\n')
# }
# if (ncol(loadings) > 1) {
	
	# b <- as.matrix(loadings)
	# n <- nrow(loadings)
	# nf <- ncol(loadings)
	# hjsq <- diag(loadings %*% t(loadings))   # communalities
	# hj <- sqrt(hjsq)
	# bh <- loadings / (hj %*% matrix(1,1,nf))
	
	# Vtemp <- n * sum(bh^4) - sum(colSums(bh^2)^2)
	# V0 <- Vtemp
	
	# for (it in 1:30) { # Never seems to need very many iterations
	# for (i in 1:(nf-1)) { # Program cycles through 2 factors
	  # jl <- i+1    # at a time
	  # for (j in jl:nf) {
	      # xj <- loadings[,i]/hj   # notation here closely
	      # yj <- loadings[,j]/hj   # follows harman
	      # uj <- xj*xj-yj*yj
	      # vj <- 2*xj*yj
	      # A <- sum(uj)
	      # B <- sum(vj)
	      # C <- t(uj) %*% uj - t(vj) %*% vj
	      # D <- 2*t(uj)%*%vj
	      # num <- D-2*A*B/n
	      # den <- C-(A^2-B^2)/n
	      # tan4p <- num/den
	      # phi <- atan2(num,den)/4
	      # angle <- phi*180/pi
	      # if (abs(phi)>.00001) {
	          # # Xj <- cos(phi)*xj+sin(phi)*yj
	          # # Yj <- -sin(phi)*xj+cos(phi)*yj
	
	
	          # Xj <- c(cos(phi))*xj+c(sin(phi))*yj
	
	          # Yj <- c(-sin(phi))*xj+c(cos(phi))*yj
	
	
	          # bj1 <- Xj*hj
	          # bj2 <- Yj*hj
	          # b[,i] <- bj1
	          # b[,j] <- bj2
	          # loadings[,i] <- b[,i]
	          # loadings[,j] <- b[,j]
	# }}}
	
	# loadings <- b
	# bh <- loadings / (hj %*% matrix(1,1,nf))
	# Vtemp <- n * sum(bh^4) - sum(colSums(bh^2)^2)
	# V <- Vtemp
	# if (abs(V-V0) < .0001) {break
	# } else {V0 <- V}
	# }
	
	# colnames(loadings) <-  c(paste('Factor ', 1:ncol(loadings), sep=''))
	
	# if (verbose == TRUE) {
		# message('\n\nVarimax Rotated Loadings:\n')
		# print(round(loadings,2))
	# }
# }
# varimaxOutput <-  list(loadings=loadings)  
# return(invisible(varimaxOutput))
# }
