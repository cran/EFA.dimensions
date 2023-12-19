




REVERSE_CODE <- function(item, max_value = NULL) {
	
	if (is.null(max_value))  max_value <- max(item)

	item_recoded <- (min(item) + max(item)) - item

	# item_recoded <- max_value + 1 - item

	return(item_recoded)	
	
}





Cronbach.alpha <- function(data) {
	
	data <- MISSING_DROP(data)

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





omega_total <- function(data, extraction = 'ML') {

	data <- MISSING_DROP(data)

	Ncases <- nrow(data)
	
	# 2017 McNeish - Thanks Coefficient Alpha, Well Take It From Here?  formula 2, p 417

	# mlfa <- MAXLIKE_FA(data, Nfactors = 1, rotation='none', verbose=FALSE)
	
	cormat <- cor(data)

	if (extraction == 'ML') efa_output <- MAXLIKE_FA(cormat, Nfactors = 1, Ncases=Ncases)

	if (extraction == 'PAF') efa_output <- PA_FA(cormat, Nfactors = 1)
	
	loadings <- efa_output$loadingsNOROT
	
	errors <- 1 - efa_output$communalities[,2]
	
	omega_t <- sum(loadings)**2 / (sum(loadings)**2 + sum(errors))

    # RMSR
    cormat_reproduced <- loadings %*% t(loadings); diag(cormat_reproduced) <- 1
    residuals <- cormat - cormat_reproduced 
    residuals.upper <- as.matrix(residuals[upper.tri(residuals, diag = FALSE)])
    rmsr <- sqrt(mean(residuals.upper^2)) # rmr is perhaps the more common term for this stat
	
	output <- list(omega_t=omega_t, rmsr=rmsr)
	
	return(invisible(output))
}	




INTERNAL.CONSISTENCY <- function(data, extraction = 'ML', reverse_these = NULL, auto_reverse = TRUE, verbose=TRUE, factormodel) {

  # deprecated  
  if (!missing(factormodel))  extraction <- factormodel

	data <- MISSING_DROP(data)

	cnoms <- colnames(data) 

	Nitems <- ncol(data)

	new_data <- data
	
	
	# reverse code these items
	if (!is.null(reverse_these)) {
		
		for (lupe in 1:length(reverse_these)) {
			
			new_data[,reverse_these[lupe]] <- REVERSE_CODE(item = data[,reverse_these[lupe]], max_value = NULL)						
				
			message('\nItem ', cnoms[lupe], ' has been reverse-coded, as requested')
			
			colnames(new_data)[colnames(new_data) == reverse_these[lupe]] <- paste(reverse_these[lupe], "rev", sep="_")
		}
	}
	

	# reverse code every item that has a negative loading on the first principal component, BUT NOT IF 
	# all of the loadings are negative 
	if (auto_reverse & is.null(reverse_these)) {
				
		pc1 <- PCA(data, Nfactors = 1, rotation = 'none', verbose = FALSE)$loadingsNOROT

		if (!all(pc1 < 0) | !all(pc1 < 0) ) {
 
	 		for (lupe in 1:ncol(data))  {
				
				if (pc1[lupe,1] < 0)  {
					
					new_data[,lupe] <- REVERSE_CODE(item = data[,lupe], max_value = NULL)	
					
					message('\nItem ', cnoms[lupe], 
					        ' has been reverse-coded due to a negative loading on the first principal component')
					        
					colnames(new_data)[lupe] <- paste(colnames(new_data)[lupe], "rev", sep="_")								
				}
			}
		}
	}

	item_noms <- colnames(new_data)
		
	if (Nitems < 3) 			
		int.consist_scale <- cbind(NA, Cronbach.alpha(new_data), NA)

	if (Nitems > 2) {
		omega_res <- omega_total(new_data, extraction=extraction)
		
	 	int.consist_scale <- cbind(omega_res$omega_t, Cronbach.alpha(new_data), omega_res$rmsr)
	}

	scale_tot <- rowSums(new_data)
	
	int.consist_dropped <- matrix(NA, Nitems, 7)
	item_stats <- matrix(NA, Nitems, 4)

	item_values <- sort(unique(as.vector(as.matrix(new_data))), decreasing=FALSE)

	resp_opt_props <- data.frame( matrix(NA, 1, length(item_values))); colnames(resp_opt_props) <- item_values	
	resp_opt_freqs <- resp_opt_props

	for (lupe in 1:Nitems) {

		if (Nitems == 2) 			
			int.consist_dropped[lupe,] <- cbind(NA, NA, NA, NA, NA, NA, NA)

		if (Nitems == 3) {		
			r_corxtd_item_total <- cor(rowSums(new_data[,-lupe]), new_data[,lupe] )  

			int.consist_dropped[lupe,] <- cbind(NA, Cronbach.alpha(new_data[,-lupe]), NA, r_corxtd_item_total)

	# colnames(int.consist_dropped) <- c('omega','alpha','alpha.z','r_mean','r_median','rmsr','r_corxed_item_total')
		}


		if (Nitems > 3) {
			omega_res <- omega_total(new_data[,-lupe], extraction=extraction)
			
			r_corxtd_item_total <- cor(rowSums(new_data[,-lupe]), new_data[,lupe] )  
	
			int.consist_dropped[lupe,] <- cbind(omega_res$omega_t, Cronbach.alpha(new_data[,-lupe]), omega_res$rmsr, r_corxtd_item_total)
		}

		item_dat <- new_data[,lupe]

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

		freqs[,paste(missingcols,  sep="")] = rep(0, length(missingcols))
        resp_opt_freqs <- rbind.data.frame(resp_opt_freqs, freqs)
		
		resp_opt_props_temp <- freqs / sum(freqs, na.rm=TRUE)
		# names(resp_opt_props_temp) <- c(0,3,2,4,1);  resp_opt_props_temp	# to verify rbinding by colname	
		resp_opt_props <- rbind.data.frame(resp_opt_props, resp_opt_props_temp)		
	}
	
	dimnames(int.consist_scale) <- list(rep("", dim(int.consist_scale)[1]))
	colnames(int.consist_scale)   <- c('omega','alpha','alpha.z','r_mean','r_median','rmsr')
	colnames(int.consist_dropped) <- c('omega','alpha','alpha.z','r_mean','r_median','rmsr','r_corxed_item_total')

	colnames(item_stats) <- c('N','Mean','SD','item_total_r')

    resp_opt_freqs <- resp_opt_freqs[-1,] 
    resp_opt_props <- resp_opt_props[-1,] 

	rownames(resp_opt_freqs) <- rownames(resp_opt_props) <- rownames(item_stats) <- rownames(int.consist_dropped) <- item_noms


	if (verbose) {

		if (Nitems < 4) 			
		message('\nNA values below indicate that a statistic could not be computed.')	
		
		message('\n\nReliability, interitem correlations, & 1-factor model fit (rmsr):\n')	
		print(round(int.consist_scale,2), print.gap=4)
		
		message('\n\nReliability, correlations, & 1-factor model fit (rmsr) when an item is dropped:\n')	
		print(round(int.consist_dropped,2), print.gap=4)
		
		message('\n\nItem statistics:\n')	
		print(round(item_stats,2), print.gap=4)
		
		message('\n\nResponse option frequencies:\n')	
		print(resp_opt_freqs, print.gap=4)
		
		message('\n\nResponse option proportions:\n')	
		print(round(resp_opt_props,2), print.gap=4)
			
	}

	output <- list(int.consist_scale=int.consist_scale, int.consist_dropped=int.consist_dropped, item_stats=item_stats,
	               resp_opt_freqs=resp_opt_freqs, resp_opt_props=resp_opt_props, new_data=new_data)
	
	return(invisible(output))

}



