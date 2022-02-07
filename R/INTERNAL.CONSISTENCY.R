

Cronbach.alpha <- function(data) {

	# Reliability -- see esp the Footnote, at the bottom.pdf
	itemSDs <- apply(data, 2, sd)

	totSD <- sd(rowSums(data))

	Nitems <- ncol(data)

	Calpha <- (Nitems / (Nitems - 1)) * (1 - ( sum(itemSDs**2) / totSD**2))

	cormat <- cor(data)
	
	Calpha.z <- (1 - Nitems / sum(cormat)) * (Nitems / (Nitems - 1))

	r_mean <- mean(cormat[lower.tri(cormat)])
	
	r_median <- median(cormat[lower.tri(cormat)])
	
	output <- cbind(Calpha, Calpha.z, r_mean, r_median)
	
	return(invisible(output))
}





omega_total <- function(data, factormodel = 'ML') {

	# 2017 McNeish - Thanks Coefficient Alpha, Well Take It From Here?  formula 2, p 417

	# mlfa <- MAXLIKE_FA(data, Nfactors = 1, rotate='none', verbose=FALSE)

	if (factormodel == 'ML') efa_output <- MAXLIKE_FA(data, Nfactors = 1, rotate='none', verbose=FALSE)

	if (factormodel == 'PAF') efa_output <- PA_FA(data, Nfactors = 1, rotate='none', verbose=FALSE)
	
	loadings <- efa_output$loadingsNOROT
	
	errors <- efa_output$uniquenesses
	
	omega_t <- sum(loadings)**2 / (sum(loadings)**2 + sum(errors))

	rmsr <- efa_output$fit_coefficients$RMSR
	
	output <- list(omega_t=omega_t, rmsr=rmsr)
	
	return(invisible(output))
}	






INTERNAL.CONSISTENCY <- function(data, factormodel = 'ML', verbose=TRUE) {

	Nitems <- ncol(data)

	item_noms <- colnames(data)

	omega_res <- omega_total(data, factormodel=factormodel)
		
	int.consist_scale <- cbind(omega_res$omega_t, Cronbach.alpha(data), omega_res$rmsr)

	scale_tot <- rowSums(data)
	
	int.consist_dropped    <- matrix(NA, Nitems, 6)
	item_stats <- matrix(NA, Nitems, 4)
	freqs <- matrix(NA, Nitems, )

	item_values <- sort(unique(as.vector(as.matrix(data))), decreasing=FALSE)

	resp_opt_props <- data.frame( matrix(NA, 1, length(item_values))); colnames(resp_opt_props) <- item_values	
	resp_opt_freqs <- resp_opt_props

	for (lupe in 1:Nitems) {

		omega_res <-omega_total(data[,-lupe])
		
		int.consist_dropped[lupe,] <- cbind(omega_res$omega_t, Cronbach.alpha(data[,-lupe]), omega_res$rmsr)  

		item_dat <- data[,lupe]

		item_N  <- length(na.omit(item_dat))

		item_MN <- mean(na.omit(item_dat))

		item_SD <- sd(na.omit(item_dat))
		
		item_r <- cor(cbind(item_dat, scale_tot))[1,2]
				
		item_stats[lupe,] <- cbind(item_N, item_MN, item_SD, item_r)
				
		freqs <- table(item_dat)
		noms <- names(freqs)
		freqs <- data.frame(matrix(unlist(freqs), nrow=1))  
		colnames(freqs) <- noms

		# if freq = 0 for a response option, create a column for it & enter NA
		missingcols <- setdiff(colnames(resp_opt_freqs), colnames(freqs))

		freqs[,paste(missingcols,  sep="")] = rep(NA, length(missingcols))
        resp_opt_freqs <- rbind.data.frame(resp_opt_freqs, freqs)
		
		resp_opt_props_temp <- freqs / sum(freqs, na.rm=TRUE)
		# names(resp_opt_props_temp) <- c(0,3,2,4,1);  resp_opt_props_temp	# to verify rbinding by colname	
		resp_opt_props <- rbind.data.frame(resp_opt_props, resp_opt_props_temp)
		
	}
	
	dimnames(int.consist_scale) <- list(rep("", dim(int.consist_scale)[1]))
	colnames(int.consist_scale) <- colnames(int.consist_dropped) <- c('omega','alpha','alpha.z','r_mean','r_median','rmsr')

	colnames(item_stats) <- c('N','Mean','SD','item_total_r')

    resp_opt_freqs <- resp_opt_freqs[-1,] 
    resp_opt_props <- resp_opt_props[-1,] 

	rownames(resp_opt_freqs) <- rownames(resp_opt_props) <- rownames(item_stats) <- rownames(int.consist_dropped) <- item_noms


	if (verbose) {
		
		message('\n\nReliability, interitem correlations, & 1-factor model fit (rmsr):\n')	
		print(round(int.consist_scale,2), print.gap=4)
		
		message('\n\nReliability, interitem correlations, & 1-factor model fit (rmsr) when an item is int.consist_dropped:\n')	
		print(round(int.consist_dropped,2), print.gap=4)
		
		message('\n\nItem statistics:\n')	
		print(round(item_stats,2), print.gap=4)
		
		message('\n\nResponse option frequencies:\n')	
		print(resp_opt_freqs, print.gap=4)
		
		message('\n\nResponse option proportions:\n')	
		print(round(resp_opt_props,2), print.gap=4)
			
	}

	output <- list(int.consist_scale=int.consist_scale, int.consist_dropped=int.consist_dropped, item_stats=item_stats,
	               resp_opt_freqs=resp_opt_freqs, resp_opt_props=resp_opt_props)
	
	return(invisible(output))

}



