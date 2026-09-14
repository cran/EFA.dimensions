

OMEGA <- function(data, corkind = 'pearson', 
                  Nfactors = 4,
                  bifactor_kind = c('SL', 'SLiD', 'DSL', 'bifactorT', 'bigeominT'),
                  EFA_options = list(extraction = 'minres', rotation = 'oblimin'),
                  LV_options = list(group_keys = NULL, 
                                    estimator = 'MLR', 
                                    rotation = 'bigeomin',
                                    resid_correls = NULL, 
                                    LV_names = NULL,
                                    ordered = FALSE),
                  schmid_options = list(extraction = 'minres', rotation = 'oblimin'),
                  GPA_options = list(delta = .01, 
                                     epsilon = .00001, 
                                     normalize = FALSE, 
                                     maxit = 1000, 
                                     randomStarts = 50),
                  min_loading = .2, display = 2) {
  
  # auto_reverse?
  
  if (display > 0)  cat('\n\nOMEGA:') 
    
  #############################  argument checks   ##################################
  
  # is the corkind method valid?
  corkind <- corkind_check(corkind)
  
  # are the bifactor_kind valid?
  if (!all(bifactor_kind %in% c('SL', 'SLiD', 'DSL', 
                                'bifactorQ', 'bifactorT', 'bigeominQ', 'bigeominT',
                                'CFA', 'ESEM', 'none'))) {
    cat('\nThe entries for bifactor_kind contain values other than')
    cat('\nSL, SLiD, DSL, bifactorQ, bifactorT, bigeominQ, bigeominT, CFA, ESEM, none.')
    cat('\n"bifactor_kind" will be set to the default.')
    bifactor_kind <- c('SL', 'SLiD', 'DSL', 'bifactorT', 'bigeominT')
  }
  
  # check EFA_options
  EFA_options <- EFA_options_check(EFA_options)
  
  # if GPArotation will be used, EFA_options$rotation must be 'none'
  # none of the other bifactor methods use EFA
  if (any(bifactor_kind %in% c('bifactorQ', 'bifactorT', 'bigeominT', 'bigeominQ'))) {
    if (EFA_options$rotation != 'none') {
      cat('\n\nA GPArotation (one of bifactorQ, bifactorT, bigeominT, or bigeominQ) was')
      cat('\nwas specified in bifactor_kind, which means that')
      cat('\n the EFA rotation must be changed to, "none"\n') 
      EFA_options$rotation <- 'none'
    }
  }
  
  # check schmid_options if schmid will be used
  if (any(bifactor_kind %in% c('SL', 'SLiD', 'DSL'))) 
    schmid_options <- schmid_options_check(schmid_options)
  
  # if schmid will be used, is the EFA_options rotation method valid?
  if (any(bifactor_kind %in% c('SL', 'SLiD', 'DSL'))) {
    if (!EFA_options$rotation %in% c('bentlerQ', 'bentlerT', 
                                     'entropy', 'equamax', 'geominQ', 'geominT', 
                                     'oblimax', 'oblimin', 'promax', 'quartimax', 'quartimin', 
                                     'simplimax', 'varimax', 'none')) {
      cat('\nThe EFA_options entry for rotation, ', EFA_options$rotation, 
          ', is not one of the options when bifactor_kind = ', bifactor_kind, sep='')
      cat('\n"promax" will be used instead.')
      EFA_options$rotation <- 'promax'
    }
  }
  
  # check LV_options if bifactor_kind = CFA or ESEM
  if (any(bifactor_kind %in% c('CFA','ESEM')))  LV_options <- LV_options_check(LV_options)
  
  # LV_options$group_keys is required when bifactor_kind = CFA
  if (any(bifactor_kind == 'CFA' & is.null(LV_options$group_keys))) {
    cat('\n\nCFA was specified in bifactor_kind but group_keys in LV_options = NULL')
    cat('\nEither provide LV_options$group_keys or do not request CFA.')
    cat('\nbifactor_kind was changed to bifactorT, to prevent errors')
    bifactor_kind <- 'bifactorT'
  }
  
  # check GPA_options if GPA functions will be used
  if (any(bifactor_kind %in% c('bifactorQ', 'bifactorT', 'bigeominT', 'bigeominQ')))
    GPA_options <- GPA_options_check(GPA_options)
  
  
  ############################  rawdata, cormat setup   ##########################
  
  data <- MISSING_DROP(data)
  
  # # reverse code every item that has a negative loading on the first principal component, BUT NOT IF 
  # # all of the loadings are negative 
  # if (auto_reverse) {
  #   new_data <- data
  #   cnoms <- colnames(data) 
  #   
  #   pc1 <- PCA(data, Nfactors = 1, rotation = 'none', verbose = FALSE)$loadingsNOROT
  #   
  #   if (!all(pc1 < 0) | !all(pc1 < 0) ) {
  #     
  #     for (lupe in 1:ncol(data))  {
  #       
  #       if (pc1[lupe,1] < 0)  {
  #         
  #         new_data[,lupe] <- REVERSE_CODE(item = data[,lupe], max_value = NULL)	
  #         
  #         message('\nItem ', cnoms[lupe], 
  #                 ' has been reverse-coded due to a negative loading on the first principal component')
  #         
  #         colnames(new_data)[lupe] <- paste(colnames(new_data)[lupe], "rev", sep="_")								
  #       }
  #     }
  #   }
  #   data <- new_data
  # }
  
  # set up cormat
  cordat <- setupcormat(data=data, corkind=corkind, Ncases=NULL)
  cormat <- cordat$cormat
  Ncases <- cordat$Ncases
  
  
  ##################################  Nfactors   ####################################
  
  # Nfactors should be the # of group factors + 1 (for the general factor)
  
  # N_group_factors is computed & used inside bifactor_engine, not here
  
  # Nfactors should be least 2, but 3 is better
  if (Nfactors == 1) {
    cat('\nThe analyses cannot be conducted when Nfactors = 1. It will be set to 3 instead.')
    Nfactors <- 3
  }
  
  
  ###################################################################################
  
  omega_total_McD <- 
    omega_total_SL <-           omega_hierl_SL <- 
    omega_total_SLiD <-         omega_hierl_SLiD <- 
    omega_total_DSL <-          omega_hierl_DSL <-           
    omega_total_bifactorT <-    omega_hierl_bifactorT <-      
    omega_total_bifactorQ <-    omega_hierl_bifactorQ <-      
    omega_total_bigeominT <-    omega_hierl_bigeominT <-      
    omega_total_bigeominQ <-    omega_hierl_bigeominQ <- 
    omega_total_ESEM <-         omega_hierl_ESEM <- 
    loadings_SL <-     
    loadings_SLiD <- 
    loadings_DSL <-    
    loadings_bifactorT <-  
    loadings_bifactorQ <-  
    loadings_bigeominT <-  
    loadings_bigeominQ <-  
    loadings_ESEM <-  
    loadingsNOROT <- NULL
  
  outpmat <- c()
  
  
  # always provide McDonald's omega
  # McDonald's omega - McNeish p 417 formula 2  -- using 1-factor EFA & no bifactor/S-L
  efa_output <- EFA.dimensions::EFA(data=cormat, 
                                    extraction = EFA_options$extraction, 
                                    rotation='none', 
                                    corkind=corkind, Ncases=Ncases, Nfactors = 1, 
                                    verbose=FALSE)
  loadings_McD <- efa_output$loadingsNOROT
  errors <- 1 - efa_output$communalities
  omega_total_McD <- sum(loadings_McD)**2 / (sum(loadings_McD)**2 + sum(errors))
  rmsr <- RMSR_boc(cormat, cormat_reproduced = reproduced_R(loadings_McD))
  outpmat <- rbind(outpmat, cbind(omega_total_McD, NA, NA, NA, NA, NA, rmsr, NA))
  
  
  bif_outp_SL <- bif_outp_SLiD <- bif_outp_DSL <- bif_outp_bifactorQ <-
    bif_outp_bifactorT <- bif_outp_bigeominQ <- bif_outp_bigeominT <- bif_outp_ESEM <- NULL
  
  if (!'none' %in% bifactor_kind) {
    
    for (lupe in 1:length(bifactor_kind)) {
      
      # need unrotated loadings for 'bifactorT', 'bifactorQ', 'bigeominT', 'bigeominQ'
      if (bifactor_kind[lupe] %in% c('bifactorT', 'bifactorQ', 'bigeominT', 'bigeominQ')) {
        
        loadingsNOROT <- EFA(data=cormat, 
                             Nfactors = Nfactors, 
                             extraction = EFA_options$extraction, 
                             rotation = 'none', 
                             corkind = corkind, Ncases = Ncases, 
                             iterpaf=100, ppower = 3, 
                             verbose=FALSE)$loadingsNOROT 
      }
      
      bif_outp <- bifactor_engine(loadings = loadingsNOROT, 
                                  cormat = cormat, 
                                  rawdata = data,
                                  corkind = corkind, Ncases = Ncases, 
                                  Nfactors = Nfactors,
                                  bifactor_kind = bifactor_kind[lupe],
                                  LV_options = LV_options,
                                  schmid_options = schmid_options,
                                  min_loading = min_loading)  

      assign(paste("bif_outp_", bifactor_kind[lupe], sep=""),    bif_outp)
      
      assign(paste("loadings_", bifactor_kind[lupe], sep=""),    bif_outp$loadings_BIF)
      
      assign(paste("omega_total_", bifactor_kind[lupe], sep=""), bif_outp$omega_total)
      
      assign(paste("omega_hierl_", bifactor_kind[lupe], sep=""), bif_outp$omega_hierl)
      
      fitcoefs <- c(bif_outp$ECV_SS[1], bif_outp$coef_H[1], bif_outp$FD[1], 
                    bif_outp$ARPB$ARPB_total, bif_outp$rmsr, bif_outp$rmsr_gen)
      
      outpmat <- rbind(outpmat, cbind(bif_outp$omega_total[1], bif_outp$omega_hierl[1], t(fitcoefs)))
    }
  }
  
  
  outpmat <- as.data.frame(outpmat)
  
  if ( 'none' %in% bifactor_kind)  rownames(outpmat) <- 'McDonald'
  if (!'none' %in% bifactor_kind)  rownames(outpmat) <- c('McDonald', bifactor_kind)
  
    
  # rownames(outpmat) <- c('McD', bifactor_kind)
  colnames(outpmat) <- c('omega-T', 'omega-H', 'ECV', 'H', 'FD', 'ARPB', 'rmsr', 'rmsr-g')
  
  
  #############################  output   ###########################################

    if (display > 0) {
    
    cat('\n\nNfactors = ', Nfactors)
    cat('\nNcases = ', Ncases) 
    cat('\nextraction = ', EFA_options$extraction)
    cat('\nrotation = ', EFA_options$rotation)
    cat('\ncorkind = ', corkind)
    
    cat('\n\nOmegas & bifactor model statistics:\n\n')
    print(round(outpmat,2), print.gap=3)
  }
  
  if (display == 2) {
    
    if ('SL' %in% bifactor_kind) {
      cat('\n\n\nSchmid-Leiman loadings:\n\n')
      print(round(loadings_SL,3), print.gap=4)
      show_bifactor_stats(bif_outp_SL)
    }
    
    if ('SLiD' %in% bifactor_kind) {
      cat('\n\n\nSLiD loadings:\n\n')
      print(round(loadings_SLiD,3), print.gap=4)
      show_bifactor_stats(bif_outp_SLiD)
    }
    
    if ('DSL' %in% bifactor_kind) {
      cat('\n\n\nDirect S-L loadings:\n\n')
      print(round(loadings_DSL,3), print.gap=4)
      show_bifactor_stats(bif_outp_DSL)
    }
    
    if ('bifactorT' %in% bifactor_kind) {
      cat('\n\n\nbifactorT loadings:\n\n')
      print(round(loadings_bifactorT,3), print.gap=4)
      show_bifactor_stats(bif_outp_bifactorT)
    }
    
    if ('bifactorQ' %in% bifactor_kind) {
      cat('\n\n\nbifactorQ loadings:\n\n')
      print(round(loadings_bifactorQ,3), print.gap=4)
      show_bifactor_stats(bif_outp_bifactorQ)
    }
    
    if ('bigeominT' %in% bifactor_kind) {
      cat('\n\n\nbigeominT loadings:\n\n')
      print(round(loadings_bigeominT,3), print.gap=4)
      show_bifactor_stats(bif_outp_bigeominT)
    }
    
    if ('bigeominQ' %in% bifactor_kind) {
      cat('\n\n\nbigeominQ loadings:\n\n')
      print(round(loadings_bigeominQ,3), print.gap=4)
      show_bifactor_stats(bif_outp_bigeominQ)
    }

    if ('ESEM' %in% bifactor_kind) {
      cat('\n\n\nESEM loadings:\n\n')
      print(round(loadings_ESEM,3), print.gap=4)
      show_bifactor_stats(bif_outp_ESEM)
    }
  }
  
  output <- list(omega_total_McD       = omega_total_McD,
                 omega_total_SL        = omega_total_SL,          omega_hierl_SL         = omega_hierl_SL,
                 omega_total_SLiD      = omega_total_SLiD,        omega_hierl_SLiD       = omega_hierl_SLiD,
                 omega_total_DSL       = omega_total_DSL,         omega_hierl_DSL        = omega_hierl_DSL,           
                 omega_total_bifactorT = omega_total_bifactorT,   omega_hierl_bifactorT  = omega_hierl_bifactorT,     
                 omega_total_bifactorQ = omega_total_bifactorQ,   omega_hierl_bifactorQ  = omega_hierl_bifactorQ,     
                 omega_total_bigeominT = omega_total_bigeominT,   omega_hierl_bigeominT  = omega_hierl_bigeominT,     
                 omega_total_bigeominQ = omega_total_bigeominQ,   omega_hierl_bigeominQ  = omega_hierl_bigeominQ,
                 omega_total_ESEM      = omega_total_ESEM,        omega_hierl_ESEM       = omega_hierl_ESEM,
                 loadings_SL           = loadings_SL,    
                 loadings_SLiD         = loadings_SLiD,
                 loadings_DSL          = loadings_DSL,   
                 loadings_bifactorT    = loadings_bifactorT, 
                 loadings_bifactorQ    = loadings_bifactorQ, 
                 loadings_bigeominT    = loadings_bigeominT, 
                 loadings_bigeominQ    = loadings_bigeominQ, 
                 loadings_ESEM         = loadings_ESEM,
                 cormat = cormat,
                 Ncases = Ncases,
                 outpmat = outpmat) 
  
  return(invisible(output))
}	



