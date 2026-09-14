


PLOT_Invariance <- function(model_object, 
                            invar_model = 'Scalar',
                            plot_types = c('loadings', 'ints_slopes', 'LV_distribs'),
                            loadings_std = TRUE,
                            loadings_by_factor = FALSE) {
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(ask=TRUE)
  
  # check on invar_model
  if (!invar_model %in% c('Configural', 'Metric', 'Scalar', 'Strict', 
                         'LV_vars', 'LV_covars', 'LV_means')) {
    cat('\n\nThe entry for the invar_model argument,', invar_model, 'is not one of the options.')
    cat('\nScalar will be used instead.')
    invar_model = 'Scalar'
  }
  
  # check on plot_types
  if (!all(plot_types %in% c('loadings', 'ints_slopes', 'LV_distribs'))) {
    cat('\n\nOne or more values for the plot_types argument is not one of the options.')
    cat('\nLV_distribs will be used instead.')
    plot_types = 'LV_distribs'
  }
  
  if (inherits(model_object, 'Factorial_Invariance'))
    model_object <- eval(parse(text=paste('model_object$mod_', invar_model, sep='')))
    
  # the group variable labels/names 
  group_noms <- lavInspect(model_object, "group.label")
  
  group_num_ids <- 1:length(group_noms)
  
  # the LV names 
  LV_noms <- lavNames(model_object, type = "lv")
  
  N_LVs <- length(LV_noms)
  
  N_groups <- length(group_noms)
  
  cols <- c('blue', 'red', 'cyan2', 'darkviolet', 'chartreuse1', 'yellow',
            'burlywood3','darkseagreen1', 
            'mediumvioletred', 'darkgreen','bisque','cyan3', 'deeppink4')
  
  output <- list(loadings = NULL, estimates = NULL, LV_scores = NULL)
  
  
  if ('loadings' %in% plot_types) {
    
    # extract loadings
    if ( loadings_std) loadings <- standardizedSolution(model_object)
    if (!loadings_std) loadings <- parameterEstimates(model_object)
    
    # filter for factor loadings (op == "=~")
    loadings <-loadings[loadings$op == "=~", ]
    # print(loadings)
    
    output$loadings <- loadings
    
    grp_combins <- combn(group_num_ids, m = 2, simplify = FALSE)
    
    if (!loadings_by_factor) {
      
      for (lupe_grp_combins in 1:length(grp_combins)) {
        
        # subset and reshape loadings for Group 1 vs Group 2
        # filter for '=~' operator representing factor loadings
        loadings_g1 <- subset(loadings, loadings$op == "=~" 
                              & loadings$group == grp_combins[[lupe_grp_combins]][1])
        loadings_g2 <- subset(loadings, loadings$op == "=~" 
                              & loadings$group == grp_combins[[lupe_grp_combins]][2])
        
        g1_label <- group_noms[ grp_combins[[lupe_grp_combins]][1] ]
        g2_label <- group_noms[ grp_combins[[lupe_grp_combins]][2] ]
        titre <- paste('DIF Check for All Factor Loadings:', g1_label, ' vs. ', g2_label)
        
        plot(loadings_g1$est, loadings_g2$est, pch = 19, col = "blue", 
             xlab = paste('Group "', g1_label, '" Loadings', sep=''), 
             ylab = paste('Group "', g2_label, '" Loadings', sep=''),
             main = titre,
             xlim = c(.5, 1),
             ylim = c(.5, 1))
        
        # the diagonal line (intercept = 0, slope = 1)
        abline(a = 0, b = 1, col = "red",  lwd = 1) 
        
        # text labels (pos = 3 places text above the point)
        text(loadings_g1$est, loadings_g2$est, labels = loadings_g1$rhs, pos = 3, cex = 0.7) 
      }
    }
    
    if ( loadings_by_factor) {
      
      for (lupe_LVs in 1:N_LVs) {
        
        for (lupe_grp_combins in 1:length(grp_combins)) {
          
          # subset and reshape loadings for Group 1 vs Group 2
          # filter for '=~' operator representing factor loadings
          
          loadings_g1 <- subset(loadings, loadings$op == "=~" 
                                & loadings$group == which(group_num_ids == grp_combins[[lupe_grp_combins]][1])  
                                & loadings$lhs == LV_noms[lupe_LVs])
          
          loadings_g2 <- subset(loadings, loadings$op == "=~" 
                                & loadings$group == which(group_num_ids == grp_combins[[lupe_grp_combins]][2]) 
                                & loadings$lhs == LV_noms[lupe_LVs])
          
          
          # loadings_g1 <- subset(loadings, loadings$op == "=~" 
          #                       & loadings$group == which(group_noms == grp_combins[[lupe_grp_combins]][1])  
          #                       & loadings$lhs == LV_noms[lupe_LVs])
          # 
          # loadings_g2 <- subset(loadings, loadings$op == "=~" 
          #                       & loadings$group == which(group_noms == grp_combins[[lupe_grp_combins]][2]) 
          #                       & loadings$lhs == LV_noms[lupe_LVs])
          
          g1_label <- group_noms[ grp_combins[[lupe_grp_combins]][1] ]
          g2_label <- group_noms[ grp_combins[[lupe_grp_combins]][2] ]
          titre <- paste('DIF Check for',  LV_noms[lupe_LVs], ':', g1_label, ' vs. ', g2_label)
          
          min_lim <- min(c(loadings_g1$est, loadings_g2$est))
          min_lim <- min_lim - .1
          
          plot(loadings_g1$est, loadings_g2$est, pch = 19, col = "blue", 
               xlab = paste('Group "', g1_label, '" Loadings', sep=''), 
               ylab = paste('Group "', g2_label, '" Loadings', sep=''),
               main = titre,
               xlim = c(min_lim, 1),
               ylim = c(min_lim, 1))
          
          # the diagonal line (intercept = 0, slope = 1)
          abline(a = 0, b = 1, col = "red",  lwd = 1) 
          
          # text labels (pos = 3 places text above the point)
          text(loadings_g1$est, loadings_g2$est, labels = loadings_g1$rhs, pos = 3, cex = 0.7) 
        }
      }
    }
  }
  
  
  if ('ints_slopes' %in% plot_types)  {
    
    ests <- parameterEstimates(model_object)
    
    output$estimates <- ests
    
    grp1_ints <- subset(ests, ests$op == "~1" & ests$group == 1)
    grp1_slos <- subset(ests, ests$op == "=~" & ests$group == 1)
    
    # grp2_ints <- subset(ests, op == "~1" & group == 2)
    grp2_slos <- subset(ests, ests$op == "=~" & ests$group == 2)
    
    # warning if the slopes are the same for the groups 
    # (e.g., when loadings were constricted to be equal) = any model other than Configural
    if (all(round(grp1_slos$est,4) == round(grp2_slos$est,4))) {
      cat('\n\nWARNING: The slopes for the groups are apparently the same.')
      cat('\nThis occurs when the slopes are constrained to be equal.')
      cat('\nConsider using a "Configural" model instead.\n\n')
    }
    
    for (lupe_LVs in 1:N_LVs) {
      
      for (items in 1:nrow(grp1_slos)) {
        
        if (grp1_slos$lhs[items] == LV_noms[lupe_LVs]) {
          
          int1 <- grp1_ints$est[items]
          slo1 <- grp1_slos$est[items]
          
          don1 <- cbind( c(-2.5,2.5), 
                         c( (int1 + slo1 * -2.5),
                            (int1 + slo1 *  2.5) ))
          
          plot(don1, type = "l",
               main = paste('Predicted Scores for "', grp1_ints[items,'lhs'], '" by Group', sep=''),
               xlab = paste('Latent "', LV_noms[lupe_LVs], '" Scores', sep=''),
               ylab = paste('Predicted', grp1_ints[items,'lhs'], 'Score'), 
               col = "blue", 
               lwd = 2,
               xlim = c(-2.5, 2.5),
               # ylim = ylim,  # range(c(dens1$y, dens2$y))
          )
          
          for (lupe_grps in 2:N_groups)  {
            
            grp2_ints <- subset(ests, ests$op == "~1" & ests$group == lupe_grps)
            grp2_slos <- subset(ests, ests$op == "=~" & ests$group == lupe_grps)
            
            int2 <- grp2_ints$est[items]
            slo2 <- grp2_slos$est[items]
            
            don2 <- cbind( c(-2.5,2.5),
                           c( (int2 + slo2 * -2.5),
                              (int2 + slo2 *  2.5) ))

            lines(don2, col = cols[lupe_grps], lwd = 2)
          }
        }
        legend("topleft", legend = group_noms, col = cols[1:lupe_grps], 
               lty = 1, lwd = 2, bty = 'n')
      }
    }
  }
  
  
  if ('LV_distribs' %in% plot_types) {
    
    # extract latent variable scores
    LV_scores_list <- lavPredict(model_object)
    
    output$LV_scores <- LV_scores_list
    
    # N_groups <- length(LV_scores_list)
    # N_LVs <- ncol(LV_scores_list[[1]])
    
    mns <- matrix(NA, N_LVs, N_groups)
    
    for (lupe_LVs in 1:N_LVs) {
      
      # densities for each group (one LV at a time)
      dens <- list()
      for (lupe_grps in 1:N_groups) {
        dum <- density(LV_scores_list[[lupe_grps]][,lupe_LVs])
        dens[[lupe_grps]] <- cbind(dum$x, dum$y)
        
        # LV means
        mns[lupe_LVs,lupe_grps] <- mean(LV_scores_list[[lupe_grps]][,lupe_LVs])
      }
      
      xlim <- c(min(sapply(dens, function(x) min(x[, 1], na.rm = TRUE))),
                max(sapply(dens, function(x) max(x[, 1], na.rm = TRUE))) )
      
      ylim <- c(min(sapply(dens, function(x) min(x[, 2], na.rm = TRUE))),
                max(sapply(dens, function(x) max(x[, 2], na.rm = TRUE))) )
      
      plot(dens[[1]], type = 'l',
           main = paste('Latent Trait Distribution for "', LV_noms[lupe_LVs], '" by Group', sep=''),
           xlab = paste('Latent "', LV_noms[lupe_LVs], '" Scores', sep=''),
           ylab = "Density", 
           col = "blue", 
           lwd = 2,
           xlim = xlim,  # range(c(dens1$x, dens2$x)),
           ylim = ylim,  # range(c(dens1$y, dens2$y))
      )
      
      abline(v = mns[lupe_LVs,1], col = "blue", lty = 2) 
      # abline(v = median(dens[[1]][,2], na.rm = FALSE), col = "blue", lty = 2) 
      # print( median(dens[[1]][,2]) )
      
      for (lupe_grps in 2:N_groups)  {
        lines(dens[[lupe_grps]], col = cols[lupe_grps], lwd = 2)
        
        abline(v = mns[lupe_LVs,lupe_grps], col = cols[lupe_grps], lty = 2) 
        # abline(v = median(dens[[lupe_grps]][,2], na.rm = FALSE), col = cols[lupe_grps], lty = 2) 
        # print( median(dens[[lupe_grps]][,2]))
        
      }
      
      legend("topright", legend = group_noms, col = cols[1:lupe_grps], 
             lty = 1, lwd = 2, bty = 'n')
    } 
    
    # warning if the latent means are the same for the groups 
    # (e.g., when intercepts & loadings were not constricted to be equal)
    # need at least Scalar, & not Metric or Configural
    if (all(round(mns,4) == 0)) {
      cat('\n\nWARNING: The latent trait means for the groups are apparently the same.')
      cat('\nThis occurs in a model where the intercepts & loadings were not constrained to be equal.')
      cat('\nConsider using a "Scalar" model instead.\n')
    }
    
    cat('\n\nLatent variable means:\n\n')
    rownames(mns) <- LV_noms; colnames(mns) <- group_noms
    print(round(mns,4), print.gap=4)
  }
  

  return(invisible(output))
}


