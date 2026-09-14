

# lavaan->lav_options_checkvalues()




corkind_check <- function(corkind) {
  if (!corkind %in% c('pearson', 'kendall', 'spearman', 'gamma', 'polychoric')) {
    cat('\nThe entry for corkind, ', corkind, ', is not one of the options for this function.', sep='')
    cat('\n"pearson" will be used instead.')
    corkind <- 'pearson'
  }
  return(invisible(corkind))
}



schmid_options_check <- function(schmid_options) {
  
  # check the names of the provided schmid_options elements
  possible_eles <- c('extraction','rotation')
  eles <- names(schmid_options)
  if (!all(eles %in% possible_eles)) {
    cat('\n\nThe names of one or more elements in LV_options is not valid.')
    cat('\nThe possibilities are: extraction, and rotation.')
  }
  
  # update schmid_options with default values if any of the possible elements are missing
  not_there <- setdiff(possible_eles, eles)
  if (length(not_there) > 0) {
    
    if ('extraction' %in% not_there) schmid_options$extraction <- 'minres'
    
    if ('rotation' %in% not_there) schmid_options$rotation <- 'oblimin'
  }
  
  # is the schmid_options extraction method valid?
  # 3 of the bifactor methods use psych::schmid, which has restricted options
  # fm: the default is  minres. fm="pa" for principal axes, fm="pc" for principal 
  # components, fm = "minres" for minimum residual (OLS), pc="ml" for maximum likelihood
  if (!schmid_options$extraction %in% c('paf', 'minres', 'ml', 'pc')) {
    cat('\nThe schmid_options entry for extraction, ', schmid_options$extraction, 
        ', is not one of the options for this argument.', sep='')
    cat('\n"minres" will be used instead.')
    schmid_options$extraction <- 'minres'
  }
  
  # is the schmid_options rotation method valid?
  # rotate: the default, oblimin, produces somewhat more correlated factors than the 
  # alternative, simplimax. Other options include Promax (not Kaiser normalized) 
  # or promax (Promax with Kaiser normalization). See fa for possible oblique rotations.
  if (!schmid_options$rotation %in% c('oblimin','simplimax','Promax','promax', 'none')) { 
    cat('\nThe schmid_options entry for rotation, ', schmid_options$rotation, 
        ', is not one of the options for this function.', sep='')
    cat('\n"oblimin" will be used instead.')
    schmid_options$rotation <- 'oblimin'
  }
  return(invisible(schmid_options))
}




EFA_options_check <- function(EFA_options) {
  
  # check the names of the provided EFA_options elements
  possible_eles <- c('extraction','rotation')
  eles <- names(EFA_options)
  if (!all(eles %in% possible_eles)) {
    cat('\n\nThe names of one or more elements in LV_options is not valid.')
    cat('\nThe possibilities are: extraction and rotation')
  }
  
  # update EFA_options with default values if any of the possible elements are missing
  not_there <- setdiff(possible_eles, eles)
  if (length(not_there) > 0) {
    
    if ('extraction' %in% not_there) EFA_options$extraction <- 'minres'
    
    if ('rotation' %in% not_there) EFA_options$rotation <- 'oblimin'
  }
  
  # is the EFA_options extraction method valid?
  if (!EFA_options$extraction %in% c('paf', 'ml', 'image', 'minres', 'uls', 'ols', 'wls', 
                                     'gls', 'alpha', 'fullinfo')) {
    cat('\nThe EFA_options entry for extraction, ', EFA_options$extraction, 
        ', is not one of the options for this argument.', sep='')
    cat('\n"paf" will be used instead.')
    EFA_options$extraction <- 'paf'
  }

  # is the EFA_options rotation method valid?
  if (!EFA_options$rotation %in% c('bentlerQ', 'bentlerT', 'entropy', 'equamax', 
                                   'geominQ', 'geominT', 'oblimax', 'oblimin', 
                                   'promax', 'quartimax', 'quartimin', 
                                   'simplimax', 'varimax', 'none')) {
    cat('\nThe EFA_options entry for rotation, ', EFA_options$rotation, 
        ', is not one of the options for this argument.', sep='')
    cat('\n"oblimin" will be used instead.')
    EFA_options$rotation <- 'oblimin'
  }
  
  return(invisible(EFA_options))
}




LV_options_check <- function(LV_options) {

  # check the names & values of the provided LV_options elements
  possible_eles <- c('group_keys','estimator','rotation',
                     'resid_correls','LV_names','ordered')
  eles <- names(LV_options)
  if (!all(eles %in% possible_eles)) {
    cat('\n\nThe names of one or more elements in LV_options is not valid.')
    cat('\nThe possibilities are: group_keys, estimator, rotation,
                       resid_correls, LV_names, and ordered')
  }
  
  # update LV_options with default values if any of the possible elements are missing
  not_there <- setdiff(possible_eles, eles)
  if (length(not_there) > 0) {
    
    if ('group_keys' %in% not_there) LV_options$group_keys <- NULL
    
    if ('estimator' %in% not_there) LV_options$estimator <- 'ML'
    
    if ('rotation' %in% not_there) LV_options$rotation <- 'bigeomin'
    
    if ('resid_correls' %in% not_there) LV_options$resid_correls <- NULL
    
    if ('LV_names' %in% not_there) LV_options$LV_names <- NULL
    
    if ('ordered' %in% not_there) LV_options$ordered <- FALSE
  }
  
  # LV_options$rotation must be either 'bigeomin' or 'biquartimin' (for bifactor via lavaan::cfa)
  if (!LV_options$rotation %in% c('bigeomin','biquartimin')) {
    cat('\nLV_options$rotation must be either bigeomin or biquartimin. It will be changed to bigeomin.')
    LV_options$rotation <- 'bigeomin'
  }
  
  return(invisible(LV_options))
}




GPA_options_check <- function(GPA_options) {
  
  # check the names & values of the provided GPA_options elements
  possible_eles <- c('delta','epsilon','normalize','maxit','randomStarts')
  eles <- names(GPA_options)
  if (!all(eles %in% possible_eles)) {
    cat('\n\nThe names of one or more elements in GPA_options is not valid.')
    cat('\nThe possibilities are: delta, epsilon, normalize, maxit, and randomStarts')
  }
  
  # update GPA_options with default values if any of the possible elements are missing
  not_there <- setdiff(possible_eles, eles)
  if (length(not_there) > 0) {

    if ('delta' %in% not_there) GPA_options$delta <- .01
    
    if ('epsilon' %in% not_there) GPA_options$epsilon <- .00001
    
    if ('normalize' %in% not_there) GPA_options$normalize <- FALSE
    
    if ('maxit' %in% not_there) GPA_options$maxit <- 1000
    
    if ('randomStarts' %in% not_there) GPA_options$randomStarts <- 50
  }
  
  return(invisible(GPA_options))
}




rotation_func <- function(rotation, loadingsNOROT, GPA_options, ppower) {
  
  Nvars <- nrow(loadingsNOROT)
  
  loadingsROT <- structure <- pattern <- phi <- NULL
  
  # GPArotation - orthogonal rotations
  if (rotation %in% c('bentlerT', 'entropy', 'equamax', 'geominT', 'quartimax', 
                      'bigeominT', 'bifactorT')) {
    
    
    # while it is possible to call, e.g., 'bifactorT', that name will not work
    # for GPFRSorth -- the T must be removed
    # it is the orth portion of GPFRSorth that keeps in orthogonal
    if (rotation == 'bentlerT')   rotation <- 'bentler'
    if (rotation == 'geominT')    rotation <- 'geomin'
    if (rotation == 'bigeominT')  rotation <- 'bigeomin'
    if (rotation == 'bifactorT')  rotation <- 'bifactor'

    loadingsROT <- GPFRSorth(loadingsNOROT, method = rotation,
                             delta = GPA_options$delta,
                             epsilon = GPA_options$epsilon,
                             normalize = GPA_options$normalize,
                             maxit = GPA_options$maxit,
                             randomStarts = GPA_options$randomStarts)$loadings[1:Nvars,]
  }
  
  if (rotation == 'varimax') 
    loadingsROT <- VARIMAX(loadingsNOROT, verbose=FALSE)$loadingsV
  
  # GPArotation::parsimax(loadings)   not an exported object from 'namespace:GPArotation'
  
  # GPArotation - oblique rotations
  if (rotation %in% c('bentlerQ', 'geominQ', 'oblimin', 'oblimax', 'quartimin', 
                      'simplimax', 'bigeominQ', 'bifactorQ')) {
    
    # while it is possible to call, e.g., 'bifactorQ', that name will not work
    # for GPFRSooblq -- the Q must be removed
    # it is the oblq portion of GPFRSorth that keeps in oblique
    if (rotation == 'bentlerQ')   rotation <- 'bentler'
    if (rotation == 'geominQ')    rotation <- 'geomin'
    if (rotation == 'bigeominQ')  rotation <- 'bigeomin'
    if (rotation == 'bifactorQ')  rotation <- 'bifactor'
    
    outp <- GPFRSoblq(loadingsNOROT, method = rotation,
                      delta = GPA_options$delta,
                      epsilon = GPA_options$epsilon,
                      normalize = GPA_options$normalize,
                      maxit = GPA_options$maxit,
                      randomStarts = GPA_options$randomStarts)
    pattern <- outp$loadings[1:Nvars,]
    phi <- outp$Phi
    structure <- pattern %*% phi
  }
  
  if (rotation == 'promax' | rotation == 'PROMAX') {
    promaxOutput <- PROMAX(loadingsNOROT, ppower=ppower, verbose=FALSE)
    pattern <- promaxOutput$pattern
    structure <- promaxOutput$structure
    phi <- promaxOutput$phi
  }
  
 output <- list(loadingsROT = loadingsROT, structure =structure, pattern = pattern, phi = phi)
  
 return(invisible(output))

  # cmd_string <- paste('GPArotation::', 
  #                     rotation, 
  #                     '(A = loadings,
  #                         delta = GPA_options$delta,
  #                         epsilon = GPA_options$epsilon,
  #                         normalize = GPA_options$normalize,
  #                         maxit = GPA_options$maxit,
  #                         randomStarts = GPA_options$randomStarts)', sep='')
}



ordered_estimator_check <- function(ordered, estimator) {
  
  # make sure ordered is compatible with estimator
  
  flag <- FALSE
  
  if ( is.logical(ordered)) { if (ordered)  flag <- TRUE}
       
  if (!is.logical(ordered))  flag <- TRUE
    
  if (flag & !estimator %in% c('WLSMV','WLSM','DWLS','ULSMV','ULSM','ULS','PML')) {
      cat('\n\nThe estimator, ', estimator, ', cannot be used for ordered data.', sep='')
      cat('\nThe estimator has therefore been changed to WLSMV.\n')
      estimator <- 'WLSMV'
    }
  return(invisible(estimator))
}



ordered_data_check <- function(ordered, rawdata) {
  
  Nlevels <- sapply(rawdata, function(x) length(unique(x)))

  if (is.logical(ordered)) {
    # if ordered = TRUE
    if (ordered) {
      # if all have > 10 levels, change to ordered = FALSE
      if (all(Nlevels > 10))  {
        cat('\n\nAll variables had more than 10 levels/values\n\n')
        print(Nlevels)
        cat('\nordered was therefore changed to FALSE\n')
        ordered <- FALSE
      }
      # if some have > 10 levels, keep them as ordered
      if (ordered & any(Nlevels > 10)) {
        cat('\n\nSome variables had more than 10 levels/values & others had less than 10 levels/values.\n')
        cat('Only the variables with less than 10 levels/values were treated as ordered in the analyses.\n')
        ordered <- names(which(Nlevels <= 10))
      }
    }

    # if ordered = FALSE & any have 10 or fewer values
    if (!ordered) {
      if (any(Nlevels <= 10)) {
        cat('\n\nThese variables have 10 or fewer levels/values:\n\n')
        print(Nlevels[which(Nlevels <= 10)])
        cat('\nConsider changing ordered to TRUE in such cases.\n\n')
      }
    }
  }
  
  # if ordered = a vector of names, check the number of values for each variable 
    if (!is.logical(ordered)) {
      
      # check if the ordered variables have > 10 values
      if (any(Nlevels[ordered] > 10)) {
        cat('\n\nThese variables have more than 10 levels/values:\n\n')
        # print(names(which(Nlevels > 10)))
        print(Nlevels[which(Nlevels > 10)])
        cat('\nConsider removing their "ordered" status.\n\n')
      }
      
      # the NON ordered variables, if any
      # check if the NON ordered variables have <= 10 values
      vars_NON_ordrd <- names(rawdata)[which(!names(rawdata) %in% ordered)]
      Nlevels_NON_ordrd <- Nlevels[vars_NON_ordrd]
      if (length(Nlevels_NON_ordrd) > 0) {
        if (any(Nlevels_NON_ordrd <= 10)) {
          cat('\n\nThese variables have 10 or fewer levels/values:\n\n')
          print(Nlevels_NON_ordrd[Nlevels_NON_ordrd <= 10])
          cat('\nConsider changing their status to "ordered".\n\n')
        }
      }
      
      # vars_NON_ordrd <- names(rawdata)[which(!names(rawdata) %in% ordered)]
      # if (length(vars_NON_ordrd) > 0) {
      #   if (any(Nlevels[vars_NON_ordrd] <= 10)) {
      #     cat('\n\nThese variables have 10 or fewer levels/values:\n\n')
      #     print(Nlevels[which(Nlevels[vars_NON_ordrd] <= 10)])
      #     cat('\nConsider changing their status to "ordered".\n\n')
      #   }
      # }
    }
  return(invisible(ordered))
}



show_lavaan_stats <- function(lavaan_output, these = 'all') {
  
  if ('all' %in% these | 'fits' %in% these) {
    
    fits <- round(data.frame(lavaan::fitmeasures(lavaan_output)),2) 
    cat('\nFit coefficients\n')
    cat('\n RMSR  =',     format(fits['rmr',], nsmall = 2),   
        '  GFI =',        format(fits['gfi',], nsmall = 2),   
        '  NFI  =',       format(fits['nfi',], nsmall = 2),  
        '  BIC =',        format(fits['bic',], nsmall = 2),
        '  Chisq  =',     format(fits['chisq',], nsmall = 2),
        '\n RMSEA =',     format(fits['rmsea',], nsmall = 2), 
        '  TLI =',        format(fits['tli',], nsmall = 2),   
        '  NNFI =',       format(fits['nnfi',], nsmall = 2), 
        '  AIC =',        format(fits['aic',], nsmall = 2),
        '  df     =',     format(fits['df',], nsmall = 2),
        '\n SRMR  =',     format(fits['srmr',], nsmall = 2),  
        '  CFI =',        format(fits['cfi',], nsmall = 2),   
        '  IFI  =' ,      format(fits['ifi',], nsmall = 2), 
        '  MFI =',        format(fits['mfi',], nsmall = 2),
        '      pvalue =', format(fits['pvalue',], nsmall = 2)
    )
    if (any(grepl("robust", rownames(fits)))) {
      cat('\n\nRobust fit coefficients\n')
      cat('\n CFI robust   =', format(fits['cfi.robust',], nsmall = 2),   
          '  TLI robust =',    format(fits['tli.robust',], nsmall = 2),   
          '  NNFI robust =',  format(fits['nnfi.robust',], nsmall = 2),  
          '\n RMSEA robust =', format(fits['rmsea.robust',], nsmall = 2), 
          '  RNI robust =',    format(fits['rni.robust',], nsmall = 2),   
          '  GFI robust  =',   format(fits['gfi.robust',], nsmall = 2)
      )
    }
    
    # [1] "npar"                          "fmin"                          "chisq"                         "df"                           
    # [5] "pvalue"                        "chisq.scaled"                  "df.scaled"                     "pvalue.scaled"                
    # [9] "chisq.scaling.factor"          "baseline.chisq"                "baseline.df"                   "baseline.pvalue"              
    # [13] "baseline.chisq.scaled"         "baseline.df.scaled"            "baseline.pvalue.scaled"        "baseline.chisq.scaling.factor"
    # [17] "cfi"                           "tli"                           "cfi.scaled"                    "tli.scaled"                   
    # [21] "cfi.robust"                    "tli.robust"                    "nnfi"                          "rfi"                          
    # [25] "nfi"                           "pnfi"                          "ifi"                           "rni"                          
    # [29] "nnfi.scaled"                   "rfi.scaled"                    "nfi.scaled"                    "pnfi.scaled"                  
    # [33] "ifi.scaled"                    "rni.scaled"                    "nnfi.robust"                   "rni.robust"                   
    # [37] "logl"                          "unrestricted.logl"             "aic"                           "bic"                          
    # [41] "ntotal"                        "bic2"                          "scaling.factor.h1"             "scaling.factor.h0"            
    # [45] "rmsea"                         "rmsea.ci.lower"                "rmsea.ci.upper"                "rmsea.ci.level"               
    # [49] "rmsea.pvalue"                  "rmsea.close.h0"                "rmsea.notclose.pvalue"         "rmsea.notclose.h0"            
    # [53] "rmsea.scaled"                  "rmsea.ci.lower.scaled"         "rmsea.ci.upper.scaled"         "rmsea.pvalue.scaled"          
    # [57] "rmsea.notclose.pvalue.scaled"  "rmsea.robust"                  "rmsea.ci.lower.robust"         "rmsea.ci.upper.robust"        
    # [61] "rmsea.pvalue.robust"           "rmsea.notclose.pvalue.robust"  "rmr"                           "rmr_nomean"                   
    # [65] "srmr"                          "srmr_bentler"                  "srmr_bentler_nomean"           "crmr"                         
    # [69] "crmr_nomean"                   "srmr_mplus"                    "srmr_mplus_nomean"             "gfi"                          
    # [73] "gfi.ci.lower"                  "gfi.ci.upper"                  "gfi.ci.level"                  "gfi.robust"                   
    # [77] "gfi.ci.lower.robust"           "gfi.ci.upper.robust"           "cn_05"                         "cn_01"                        
    # [81] "gfi_lisrel"                    "agfi_lisrel"                   "pgfi"                          "mfi"                          
    # [85] "ecvi"  
  }
  
  # # unstandardized loadings
  # loadings_raw <- lavInspect(lavaan_output, what = "est")$lambda
  # cat('\nRaw loadings:\n'); print(round(loadings_raw,3))
  
  if ('all' %in% these | 'loadings_std' %in% these) {
    # standardized loadings
    loadings_std <- lavInspect(lavaan_output, what = "std")$lambda
    cat('\n\nStandardized loadings\n\n'); print(round(loadings_std,2), print.gap=4)
  }
  
  if ('all' %in% these | 'psi' %in% these) {
    # factor correlations
    psi <-  lavInspect(lavaan_output, what = "std")$psi
    cat('\nFactor correlations\n\n'); print(round(psi,2), print.gap=4)
  }
}




lavaan_model <- function(varnames, model=NULL, keys=NULL, 
                         LV_names=NULL, resid_correls=NULL,
                         Nfactors = NULL,
                         esem = FALSE,
                         bifactor = FALSE) {
  
  # for cfa & no bifactor, cfa & bifactor,   esem & no bifactor, esem & bifactor
  
  if (!is.null(model))  lav_mod <- model
  
  if (is.null(model)) {
    
    # # need keys if model is NULL
    # if (is.null(keys))
    #   stop('model and keys are both NULL. The analyses cannot proceed.')
    
    # if keys is provided, ignore the entered Nfactors & compute Nfactors from keys
    if (!is.null(keys))  Nfactors <- length(unique(keys))
    
    # if keys is NULL, then create it from the entered Nfactors
    if (is.null(keys))  keys <- 1:Nfactors

    # create a set of lavaan equations from the keys
    values <- unique(keys)
    
    # a 0 in keys indicates no group factor loading for an item
    values <- values[values != 0]
    
    if (is.null(Nfactors))  {
      stop('\nboth keys & Nfactors, in LV_options are NULL. The analyses cannot proceed.')
      # Nfactors <- length(values)
    }
    
    if (is.null(LV_names))  LV_names <- paste('LV_', 1:Nfactors, sep = '')
    
    if (!is.null(LV_names)) {
      if (length(LV_names) != Nfactors) {
        message('\n\nThe length of LV_names is not equal to the number of unique values in LV-Keys.')
        message('LV_names will be ignored.\n')
        LV_names <- NULL
      }
    }
    
    #  CFA
    if (!esem) {

      lav_mod  <- vector()
      
      if (bifactor) {
        
        dum <- paste(paste0('Gen', " =~ "), paste(varnames, collapse = ' + '))
        
        lav_mod <- paste(lav_mod, '\n', dum, collapse = '\n', sep='')
      }
      
      if (is.vector(keys)) {
        
        if (is.null(names(keys)))  names(keys) <- varnames
        
        for (lupe in 1:Nfactors) {
          
          ceci <- which(keys == values[lupe])
          
          dum <- paste(paste0(LV_names[lupe], " =~ "), 
                       paste(names(keys[ceci]), collapse = ' + '))
          
          lav_mod <- paste(lav_mod, '\n', dum, collapse = '\n', sep='')
        }
      }
      
      if (is.data.frame(keys)) {
        
        for (lupe in 1:Nfactors) {
          
          ceci <- which(keys[,2] == lupe)
          
          dum <- paste(paste0(LV_names[lupe], " =~ "), 
                       paste(keys[ceci,1], collapse = ' + '))
          
          lav_mod <- paste(lav_mod, '\n', dum, collapse = '\n', sep='')
        }
      }
      
      if (!is.null(resid_correls))
        lav_mod <- paste(lav_mod, '\n', resid_correls, collapse = '\n\n', sep='')
    }
    

    
    
    # ESEM & bifactor
    if (esem & bifactor) {
      
      lav_mod  <- paste(paste0('efa(\"BI-ESEM\")*',  'G', ' + '))
      
      for (lupe in 1:Nfactors) {
        
        if (lupe < Nfactors)
          dum <- paste(paste0('efa(\"BI-ESEM\")*',  LV_names[lupe], ' + '))
        
        if (lupe == Nfactors)
          dum <- paste(paste0('efa(\"BI-ESEM\")*',  LV_names[lupe], ' =~ \n'))
        
        lav_mod <- paste(lav_mod, '\n', dum, collapse = '\n', sep='')
      }
      
      # resid_correls ???
      
      lav_mod <- paste(lav_mod, paste(varnames, collapse = ' + '))
    }
    
    
    # ESEM (& NO bifactor)
    if (esem & !bifactor) {
      
      lav_mod  <- c()   #paste(paste0('efa(\"ESEM\")*',  'G', ' + '))
      
      for (lupe in 1:Nfactors) {
        
        if (lupe < Nfactors)
          dum <- paste(paste0('efa(\"ESEM\")*',  LV_names[lupe], ' + '))
        
        if (lupe == Nfactors)
          dum <- paste(paste0('efa(\"ESEM\")*',  LV_names[lupe], ' =~ \n'))
        
        lav_mod <- paste(lav_mod, '\n', dum, collapse = '\n', sep='')
      }
      
      # resid_correls ???
      
      lav_mod <- paste(lav_mod, paste(varnames, collapse = ' + '))
    }  

    # > fit$syntax
    # [1] "# bifactory: ESEM via lavaan native efa() block\n# Mplus equivalent: F1-F2-F3 BY x1-x9 (*1);\n\nefa(\"esem\")*F1 +\n    efa(\"esem\")*F2 +\n    efa(\"esem\")*F3 =~\n    x1 + x2 + x3 + x4 + x5 + x6 +\n    x7 + x8 + x9\n"
    # > writeLines(fit$syntax)
    # # bifactory: ESEM via lavaan native efa() block
    # # Mplus equivalent: F1-F2-F3 BY x1-x9 (*1);
    # 
    # efa("esem")*F1 +
    #   efa("esem")*F2 +
    #   efa("esem")*F3 =~
    #   x1 + x2 + x3 + x4 + x5 + x6 +
    #   x7 + x8 + x9
  }
  
  # writeLines(lav_mod)
  
    
  # test the syntax validity
  tryCatch({
    parsed_table <- lavaanify(model = lav_mod)
    # message('The lavaan model syntax is valid.')
    # message('\nLook at the parameter table structure:\n')
    # print((parsed_table[, c("lhs", "op", "rhs")]))
  }, error = function(e) {
    message("Syntax Error found:")
    print(e$message)
  })
  
  return(invisible(lav_mod))
}




# generate a target matrix from keys
target_matrix <- function(keys, varnames, bifactor) {
  
  target <- matrix(0, length(keys), max(keys))
  for (lupe in 1:length(keys))  target[lupe, keys[lupe]] <- 1
  
  if (!bifactor) colnames(target) <- paste('Factor_', 1:ncol(target), sep = '')
  
  if (bifactor) {
    target <- cbind(1, target)
    colnames(target) <- c('General', (paste('Group_', 1:(ncol(target)-1), sep = '')))
  }
  
  rownames(target) <- varnames
  
  # # check for errors in target_keys
  # if (length(target_keys) != N_obsvd_vars)
  #   cat('\nThe length of target_keys is not equal to the number of variables in data.')
  # if (min(target_keys) != 1)  cat('\nThe smallest number in target_keys is not 1')
  # if (!all(seq(min(target_keys):max(target_keys)) %in% unique((target_keys))))
  #   cat('\nNot all of the possible factor numbers appear in target_keys')
  target
}




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



# standardized alpha from a correlation matrix
Cronbach.alpha.z <- function(cormat) {
  k <- nrow(cormat)
  k / (k - 1) * (1 - sum(diag(cormat)) / sum(cormat))
}




reproduced_R <- function(loadings) {
  # reproduced correlation matrix from the factor model
  # 1. Compute communalities (sum of squared loadings per row)
  communalities <- rowSums(loadings^2)
  # 2. Compute uniquenesses
  uniquenesses <- 1 - communalities
  # 3. Create Lambda * Lambda Transpose
  Lambda_Lambda_T <- loadings %*% t(loadings)
  # 4. Create the Uniqueness diagonal matrix
  Theta <- diag(uniquenesses)
  # 5. Compute the reproduced correlation matrix
  # cormat_reproduced <- Lambda_Lambda_T + Theta
  Lambda_Lambda_T + Theta
  
  # cormat_reproduced <- loadings_SL %*% t(loadings_SL)
}



RMSR_boc <- function(cormat, cormat_reproduced) {
  residuals <- cormat - cormat_reproduced 
  residuals.upper <- as.matrix(residuals[upper.tri(residuals, diag = FALSE)])
  sqrt(mean(residuals.upper^2)) 
}




# eigenvalues for PCA, PAF, & image -- used by PARALLEL & RAWPAR
eigvals <- function(cormatrix, extraction) {
  
  if (extraction=='PCA')  evals <- eigen(cormatrix)$values
  
  if (extraction=='PAF') {
    smc <- 1 - (1 / diag(solve(cormatrix)))
    diag(cormatrix) <- smc
    evals <- eigen(cormatrix)$values 
  }
  
  if (extraction=='image') { # Gorsuch 1983, p 113; Velicer 1974, EPM, 34, 564
    d <-  diag(1 / diag(solve(cormatrix)))
    gvv <- cormatrix + d %*% solve(cormatrix) %*% d - 2 * d
    s <- sqrt(d)                  #  Velicer 1974 p 565 formula (7)
    r2 <- solve(s) %*%  gvv  %*% solve(s)  # Velicer 1974 p 565 formula (5)
    evals <- cbind(eigen(r2)$values) 
  }
  
  return(invisible(evals))
}





MISSING_DROP <- function(data) {
  
  if (anyNA(data) ) {
    
    # total # of NAs
    totNAs <- sum(is.na(data))
    
    # number of rows (cases) with an NA
    nrowsNAs <- sum(apply(data, 1, anyNA))
    
    data <- na.omit(data)
    message('\nCases with missing values were found and removed from the data matrix.')
    
    message('\nThere were ', nrowsNAs, ' cases with missing values, and ', 
            totNAs, ' missing values in total.\n')		
  }
  return(invisible(data))
}




# set up cormat
setupcormat <- function(data, corkind='pearson', Ncases=NULL) {
  
  # determine whether data is a correlation matrix   
  # there is also a helper function in EFAtools = .is_cormat
  if (nrow(data) == ncol(data)) {
    if (all(diag(data==1))) {datakind = 'correlations'}} else{ datakind = 'notcorrels'}
  if (datakind == 'correlations')  {
    cormat <- as.matrix(data)
    ctype <- 'from user'
    if (is.null(Ncases)) message('\nNcases must be provided when data is a correlation matrix.\n')
  }
  
  if (datakind == 'notcorrels') {
    
    Ncases <- nrow(data)
    
    if (anyNA(data) ) {
      data <- na.omit(data)
      message('\nCases with missing values were found and removed from the data matrix.\n')
    }
    
    if (corkind=='pearson')     {cormat <- cor(data, method='pearson');  ctype <- 'Pearson'}
    if (corkind=='kendall')     {cormat <- cor(data, method='kendall');  ctype <- 'Kendall'}
    if (corkind=='spearman')    {cormat <- cor(data, method='spearman'); ctype <- 'Spearman'} 
    
    # are the data whole numbers?
    if (all((data - round(data)) == 0)) { wholenums = 1 } else { wholenums = 0 }
    
    if (corkind=='polychoric')  {
      if (wholenums == 1)  {
        cormat <- POLYCHORIC_R(data, verbose=FALSE)
        ctype <- 'polychoric'
      }
      if (wholenums == 0)  {
        cormat <- cor(data, method='pearson')
        ctype <- 'Pearson'
        message('\nPolychoric correlations were specified but there are values in the')
        message('\nreal data matrix that are not whole numbers. Pearson correlations')
        message('\nwill be used instead.') 
      }
    }
    
    if (corkind=='gamma') {
      if (wholenums == 1)  {
        cormat <- Rgamma(data, verbose=FALSE)
        ctype <- 'Goodman-Kruskal gamma'
      }
      if (wholenums == 0)  {
        cormat <- cor(data, method='pearson')
        ctype <- 'Pearson'
        message('\nGoodman-Kruskal gamma correlations were specified but there are')
        message('\nvalues in the real data matrix that are not whole numbers.')
        message('\nPearson correlations will be used instead.') 
      }
    }
  }
  
  # smooth cormat if it is not positive definite
  eigenvalues <- eigen(cormat)$values
  if (min(eigenvalues) <= 0)  cormat <- psych::cor.smooth(cormat)
  
  Output <- list(cormat=cormat, ctype=ctype, Ncases=Ncases, datakind=datakind) 
  
  return(invisible(Output))
}





# \cr\cr {Thompson, L. A. 2007. R (and S-PLUS) Manual to Accompany Agresti's Categorical Data 
# Analysis (2002) 2nd edition. https://usermanual.wiki/Document/R2020and20SPLUS2020Manual20to20Accompany20Agrestis20Categorical20Data20Analysis.540471458#google_vignette}


# Goodman-Kruskal gamma - Laura Thompson
# https://usermanual.wiki/Document/R2020and20SPLUS2020Manual20to20Accompany20Agrestis20Categorical20Data20Analysis.540471458#google_vignette
# see p 25 of
# *** 2006 Thompson - S-PLUS (and R) Manual to Accompany Agresti's Categorical Data Analysis - Splusdiscrete2

Gamma.f <- function(x)
{
  # x is a matrix of counts.  You can use output of crosstabs or xtabs.
  n <- nrow(x)
  m <- ncol(x)
  res <- numeric((n-1)*(m-1))
  for(i in 1:(n-1)) {
    for(j in 1:(m-1)) res[j+(m-1)*(i-1)] <- x[i,j]*sum(x[(i+1):n,(j+1):m])
  }
  C <- sum(res)
  res <- numeric((n-1)*(m-1))
  iter <- 0
  for(i in 1:(n-1))
    for(j in 2:m) {
      iter <- iter+1; res[iter] <- x[i,j]*sum(x[(i+1):n,1:(j-1)])
    }
  D <- sum(res)
  gamma <- (C-D)/(C+D)
}



Rgamma  <- function (donnes, verbose=TRUE) {
  
  rgamma <- matrix(-9999,ncol(donnes),ncol(donnes))
  
  for (i in 1:(ncol(donnes)-1) ) {
    for (j in (i+1):ncol(donnes) ) {
      
      dat <- donnes[,c(i,j)]
      
      rgamma[i,j] <- rgamma[j,i] <- Gamma.f(table(dat[,1], dat[,2]))
    }
  }
  diag(rgamma) <- 1
  
  if (verbose ) {
    cat("\nGoodman-Kruskal gamma correlations (Thompson, 2006): \n\n" )
    print(round(rgamma,3))
  }
  
  return(invisible(rgamma))
}





VarianceExplained <- function(eigenvalues, loadingsNOROT=NULL, loadingsROT=NULL, phi=NULL) {
  
  outp <- data.frame(matrix(-999, length(eigenvalues), 9))
  cnoms <- c('Eig', 'prop', 'cum')
  
  # no extraction
  propvar <- eigenvalues / length(eigenvalues)
  cumvar  <- cumsum(propvar) 
  totvarexpl_1 <- cbind(eigenvalues, propvar, cumvar) 
  outp[,1:3] <- totvarexpl_1
  
  # NOROT
  if (!is.null(loadingsNOROT)) {
    if (is.null(phi))  sumsqloads <- colSums(loadingsNOROT**2)
    if (!is.null(phi)) sumsqloads <- diag(phi %*% crossprod(loadingsNOROT))
    propvar <- sumsqloads / length(eigenvalues)
    cumvar  <- cumsum(propvar) 
    totvarexpl_2 <- cbind(sumsqloads, propvar, cumvar) 
    outp[1:nrow(totvarexpl_2),4:6] <- totvarexpl_2
    # cnoms <-  c(cnoms, 'SSL_NOROT', 'prop_NOROT', 'cumul_NOROT')
    cnoms <-  c(cnoms, ' Eig', ' prop', ' cum')
  }
  
  # rotated
  if (!is.null(loadingsROT)) {
    if (is.null(phi))  sumsqloads <- colSums(loadingsROT**2)
    if (!is.null(phi)) sumsqloads <- diag(phi %*% crossprod(loadingsROT))
    propvar <- sumsqloads / length(eigenvalues)
    cumvar  <- cumsum(propvar) 
    totvarexpl_3 <- cbind(sumsqloads, propvar, cumvar) 
    outp[1:nrow(totvarexpl_3),7:9] <- totvarexpl_3
    # cnoms <-  c(cnoms, 'SSL_ROT', 'prop_ROT', 'cumul_ROT')
    cnoms <-  c(cnoms, ' Eig ', ' prop ', ' cum ')
  }
  
  # Remove Columns Where All Elements are Equal 
  # outp <- data.frame(outp)
  # outp <- Filter(function(x) length(unique(x)) > 1, outp)
  
  # remove columns that are all -999
  outp2 <- round(outp,2)
  drop_these <- which(colMeans(outp) == -999)
  if (length(drop_these) > 0)  outp2 <- outp2[-drop_these]
  
  # blank-out the -999 values
  outp2[outp2 == -999] <- paste(rep(" ", 2), collapse = "")
  
  colnames(outp2) <-colnames(outp) <-  cnoms
  # colnames(outp2) <- c('SSL','prop.','cum prop', 'SSL','prop.','cum prop')
  rownames(outp2) <- rownames(outp) <- c(paste('Factor ', 1:nrow(outp2), sep=''))

  
  # # Keep columns where NOT all elements are equal to 0
  # df_clean <- outp[ , !sapply(outp, function(x) all(x == -999))]
  
  # cbind(totvarexpl_1, totvarexpl_2, totvarexpl_3)
  # outp[ , colSums(is.na(outp)) != nrow(outp)]
  # outp2 <- as.data.frame(round(outp,2))
  
  # print(outp2, print.gap=3)
  
  # varex_outp <- list(outp=outp, outp2=outp2)
  varex_outp <- outp2
  
  return(invisible(varex_outp))
}





CAF_boc <- function(cormat, cormat_reproduced=NULL) {
  
  # from 2011 Lorenzo-Seva - The Hull method for selecting the number of common factors
  
  if (is.null(cormat_reproduced)) {bigR <- cormat
  } else {bigR <- cormat - cormat_reproduced;  diag(bigR) <- 1}
  
  # smooth bigR if it is not positive definite
  if (any(eigen(bigR, symmetric = TRUE, only.values = TRUE)$values <= 0)) 
    bigR <- suppressWarnings(psych::cor.smooth(bigR))
  
  # overall KMO
  Rinv <- solve(bigR)	
  Rpart <- cov2cor(Rinv)
  cormat_sq <- bigR^2
  Rpart_sq  <- Rpart^2
  KMOnum <- sum(cormat_sq) - sum(diag(cormat_sq))
  KMOdenom <- KMOnum + (sum(Rpart_sq) - sum(diag(Rpart_sq))) 
  KMO <- KMOnum / KMOdenom
  
  CAF <- 1 - KMO
  
  return(CAF)	
}





FIT_COEFS <- function(cormat, loadings, extraction, Ncases, verbose=TRUE) {
  
  cormat_reproduced <- reproduced_R(loadings)
    
  # model statistics, based on Revelle
  # cormat_reproduced <- loadingsNOROT %*% t(loadingsNOROT); diag(cormat_reproduced) <- 1
  # model <- cormat_reproduced
  # #model <- cor.smooth(model)  #this replaces the next few lines with a slightly cleaner approach
  # #r <- cor.smooth(r)  #this makes sure that the correlation is positive semi-definite
  m.inv.r <- try(solve(cormat_reproduced, cormat), silent=TRUE)
  Nvars <- nrow(loadings)
  Nfactors <- ncol(loadings)
  dfMODEL <- Nvars * (Nvars - 1) / 2 - Nvars * Nfactors + (Nfactors * (Nfactors - 1) / 2)
  objective <- sum(diag((m.inv.r))) - log(det(m.inv.r)) - Nvars 
  chisqMODEL <- objective * ((Ncases - 1) - (2 * Nvars + 5) / 6 - (2 * Nfactors) / 3) # from Tucker & from factanal
  if(!is.nan(chisqMODEL)) if (chisqMODEL < 0) {chisqMODEL <- 0}  
  if (dfMODEL > 0) {pvalue <- pchisq(chisqMODEL, dfMODEL, lower.tail = FALSE)} else {pvalue <- NA}
  
  # the null model
  Fnull <- sum(diag((cormat))) - log(det(cormat)) - Nvars  
  chisqNULL <-  Fnull * ((Ncases - 1) - (2 * Nvars + 5) / 6 )
  dfNULL <- Nvars * (Nvars - 1) / 2
  
  
  # RMSR
  RMSR <- RMSR_boc(cormat, cormat_reproduced)
  # residuals <- cormat - cormat_reproduced 
  # residuals.upper <- as.matrix(residuals[upper.tri(residuals, diag = FALSE)])
  # mnsqdresid <- mean(residuals.upper^2) # mean of the off-diagonal squared residuals (as in Waller's MicroFact)
  # RMSR <- sqrt(mean(residuals.upper^2)) # rmr is perhaps the more common term for this stat
  # # no srmsr computation because it requires the SDs for the variables in the matrix
  
  
  # GFI (McDonald, 1999), & was also from Waller's MicroFact: 
  # 1 - mean-squared residual / mean-squared correlation
  residuals <- cormat - cormat_reproduced
  residuals.upper <- as.matrix(residuals[upper.tri(residuals, diag = FALSE)])
  mnsqdresid <- mean(residuals.upper^2) # mean of the off-diagonal squared residuals (as in Waller's MicroFact)
  mnsqdcorrel <- mean(cormat[upper.tri(cormat, diag = FALSE)]^2)
  GFI <- 1 - (mnsqdresid / mnsqdcorrel)
  
  
  # CAF from Lorenzo-Seva, Timmerman, & Kiers (2011)
  CAF <- CAF_boc(cormat, cormat_reproduced=cormat_reproduced)
  

  RMSEA <- TLI<- CFI <- MFI <- BIC <- AIC <- CAIC <- SABIC <- NA
    
  if (!extraction %in% c('PCA', 'pca','IMAGE','image')) {  

    RMSEA <- sqrt(max(((chisqMODEL - dfMODEL) / (Ncases - 1)),0) / dfMODEL)
    
    # TLI - Tucker-Lewis index (Tucker & Lewis, 1973) = 
    # NNFI - nonnormed fit index (Bentler & Bonett, 1980)
    t1 <- chisqNULL / dfNULL - chisqMODEL / dfMODEL
    t2 <- chisqNULL / dfNULL - 1 
    TLI <- 1
    if(t1 < 0 && t2 < 0) {TLI <- 1} else {TLI <- t1/t2} # lavaan   else {TLI <- 1}  
    NNFI <- TLI
    
    CFI <- ((chisqNULL - dfNULL) - (chisqMODEL - dfMODEL)) / (chisqNULL - dfNULL)
    
    # MacDonald & Marsh (1990) MFI = an absolute fit index that does not depend  
    # on comparison with another model  (T&F, 2001, p 700)
    MFI <- exp (-.5 * ( (chisqMODEL - dfMODEL) / Ncases))
    
    BIC <- chisqMODEL - dfMODEL * log(Ncases)
    
    # AIC Akaike Information Criteria (T&F, 2001, p 700)
    # not on a 0-1 scale; & the value/formula varies across software
    AIC <- chisqMODEL - 2 * dfMODEL
    
    # CAIC Consistent Akaike Information Criteria (T&F, 2001, p 700)
    # not on a 0-1 scale; & the value/formula varies across software
    CAIC <- chisqMODEL - (log(Ncases) + 1) * dfMODEL
    
    # SABIC -- Sample-Size Adjusted BIC (degree of parsimony fit index)		
    # Kenny (2020): "Like the BIC, the sample-size adjusted BIC or SABIC places a penalty 
    # for adding parameters based on sample size, but not as high a penalty as the BIC.  
    # Several recent simulation studies (Enders & Tofighi, 2008; Tofighi, & Enders, 2007) 
    SABIC <- chisqMODEL + log((Ncases+2) / 24) * (Nvars * (Nvars+1) / 2 - dfMODEL)
    
    # mirt:   SABIC <- (-2) * logLik + tmp*log((N+2)/24)
  }
  
  fitcoefsOutput <- list(cormat_reproduced = cormat_reproduced, 
                         chisqMODEL = chisqMODEL, dfMODEL = dfMODEL, pvalue = pvalue,
                         chisqNULL = chisqNULL, dfNULL = dfNULL,
                         RMSR=RMSR, GFI=GFI, CAF=CAF,
                         RMSEA=RMSEA, TLI=TLI, CFI=CFI, MFI=MFI, BIC=BIC, AIC=AIC, 
                         CAIC=CAIC, SABIC=SABIC)
  
  if (verbose) {  
    
    cat('\n\n\nFit Coefficients:')
    
    cat('\n\nChi square = ', round(chisqMODEL,2),
            '   df = ', dfMODEL,'    p = ', round(pvalue,5))
    
    cat('\n\nNull Model Chi square = ', round(chisqNULL,2), '   df = ', dfNULL)
    
    cat('\n\nRMSR = ', round(RMSR,2))
    
    cat('\n\nGFI (McDonald) = ', round(GFI,2))
    
    cat('\n\nCAF = ', round(CAF,2))
    
    if (!extraction %in% c('PCA', 'pca','IMAGE','image')) {  
      
      cat('\n\nRMSEA = ', round(RMSEA,3))
      
      cat('\n\nTLI = ', round(TLI,2))
      
      cat('\n\nCFI = ', round(CFI,2))
      
      cat('\n\nMFI = ', round(MFI,2))
      
      cat('\n\nAIC = ', round(AIC,2))
      
      cat('\n\nCAIC = ', round(CAIC,2))
      
      cat('\n\nBIC = ', round(BIC,2))
      
      cat('\n\nSABIC = ', round(SABIC,2))
    }
  }
  
  return(invisible(fitcoefsOutput))    
}     






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
    
    
    if (verbose ) {
      
      cat('\n\n\nVarimax Rotated Loadings:\n\n')
      print(round(loadingsV,2), print.gap=3)
      
      cat('\n\nThe rotation matrix:\n\n')
      print(round(rotmatV,2), print.gap=)
    }
  }
  
  varimaxOutput <-  list(loadingsNOROT=loadings, loadingsV=loadingsV, rotmatV=rotmatV, 
                         cormat_reproduced=cormat_reproduced)  
  
  return(invisible(varimaxOutput))
  
}





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
    if (verbose ) message('\nWARNING: There was only one factor. Rotation was not performed.\n')
  }
  
  if (ncol(loadings) > 1) {
    
    # varimax
    vmaxres <- varimax(loadings, normalize=TRUE)   # SPSS normalizes them
    loadingsV <- vmaxres$loadings[]	
    rotmatV <- vmaxres$rotmat
    
    promaxres <- promax(loadingsV, m=ppower)
    
    bigA <- rotmatV %*% promaxres$rotmat
    
    phi  <- solve(t(bigA) %*% bigA)
    colnames(phi) <- rownames(phi) <- c(paste('Factor ', 1:ncol(loadingsV), sep=''))
    
    Pstructure <- promaxres$loadings %*% phi  # promax structure
    Ppattern   <- promaxres$loadings[]  # promax loadings/pattern
    
    # reproduced correlation matrix
    cormat_reproduced <- Pstructure %*% t(Ppattern); diag(cormat_reproduced) <- 1
    
    
    if (verbose ) {
      
      # message('\nUnrotated Loadings:\n')
      # print(round(B,2))
      
      cat('\n\nPromax Rotation Pattern Matrix:\n\n')
      print(round(Ppattern,2), print.gap=3)
      
      cat('\n\nPromax Rotation Structure Matrix:\n\n')
      print(round(Pstructure,2), print.gap=3)
      
      cat('\n\nPromax Rotation Factor Correlations:\n\n')
      print(round(phi,2), print.gap=3)
    }
    
    promaxOutput <- list(loadingsNOROT=loadings, pattern=Ppattern, structure=Pstructure, 
                         phi=phi, cormat_reproduced=cormat_reproduced)
    
  }
  return(invisible(promaxOutput))
}






# Image Factor Extraction (Gorsuch 1983, p 113; Velicer 1974, EPM, 34, 564)

IMAGE_FA <- function (cormat, Nfactors, Ncases) {
  
  smcINITIAL <- 1 - (1 / diag(solve(cormat)))  # initial communalities
  
  eigenvalues <- eigen(cormat)$values
  
  cnoms <- colnames(cormat)
  
  # factor pattern for image analysis Velicer 1974 p 565 formula (2)
  d <-  diag(1 / diag(solve(cormat)))
  gvv <- cormat + d %*% solve(cormat) %*% d - 2 * d
  s <- sqrt(d)                     #  Velicer 1974 p 565 formula (7)
  r2 <- solve(s) %*%  gvv  %*% solve(s)    #  Velicer 1974 p 565 formula (5)
  eigval <- diag(eigen(r2) $values)
  eigvect <- eigen(r2) $vectors
  l <- eigvect[,1:Nfactors]
  dd <- sqrt(eigval[1:Nfactors,1:Nfactors])
  
  loadingsNOROT <- as.matrix(s %*% l %*% dd)      #  Velicer 1974 p 565 formula (2)
  
  communalities <- as.matrix(diag(loadingsNOROT %*% t(loadingsNOROT))) 	
  communalities <- cbind(smcINITIAL, communalities) 
  rownames(communalities) <- cnoms
  colnames(communalities) <- c('Initial', 'Extraction')
  
  imageOutput <- list(loadingsNOROT=loadingsNOROT, communalities=communalities)
  
  return(invisible(imageOutput))
}





# Maximum likelihood factor analysis - using factanal from stats

MAXLIKE_FA <- function (cormat, Nfactors, Ncases) {
  
  smcINITIAL <- 1 - (1 / diag(solve(cormat)))  # initial communalities
  
  eigenvalues <- eigen(cormat)$values
  
  Nvars <- dim(cormat)[2]
  
  cnoms <- colnames(cormat)
  
  # factanal often generates errors
  # the code below uses fa from psych when factanal produces an error
  essaye1 <- try(factanalOutput <- 
                   factanal(covmat = as.matrix(cormat), n.obs = Ncases, factors = Nfactors, 
                            rotation = 'none'), silent=TRUE)
  
  loadingsNOROT <- communalities <- NA	
  
  if (!inherits(essaye1, "try-error")) {
    
    # chisqMODEL <- unname(factanalOutput$STATISTIC)
    
    # dfMODEL <- unname(factanalOutput$dof)
    
    # pvalue <- unname(factanalOutput$PVAL)
    
    loadingsNOROT <- factanalOutput$loadings[1:dim(factanalOutput$loadings)[1],
                                             1:dim(factanalOutput$loadings)[2], drop=FALSE]
    
    # uniquenesses <- factanalOutput$uniquenesses
  }	
  
  if (inherits(essaye1, "try-error")) {
    
    # using fa from psych if factanal produces an error
    essaye2 <- try(faOutput <- fa(cormat, Nfactors, rotate="Promax", fm="mle"), silent=TRUE) 
    
    if (!inherits(essaye2, "try-error")) {	
      
      # chisqMODEL <- (Ncases - 1 - (2 * Nvars  +  5) / 6 - (2 * Nfactors) / 3) * faOutput$criteria[1] # from psych ?fa page
      
      # dfMODEL <- faOutput$dof
      
      # if (dfMODEL > 0) {pvalue <- pchisq(chisqMODEL, dfMODEL, lower.tail = FALSE)} else {pvalue <- NA}
      
      loadingsNOROT <- faOutput$Structure[1:dim(faOutput$Structure)[1],
                                          1:dim(faOutput$Structure)[2], drop=FALSE]
      
      # uniquenesses <- faOutput$uniquenesses
    }
    if (inherits(essaye2, "try-error")) 	
      message('\nerrors are produced when Nfactors = ', Nfactors, '\n')	
  }
  
  # # the null model
  # Fnull <- sum(diag((cormat))) - log(det(cormat)) - Nvars  
  
  # chisqNULL <-  Fnull * ((Ncases - 1) - (2 * Nvars + 5) / 6 )
  
  # dfNULL <- Nvars * (Nvars - 1) / 2
  
  # if there are no errors
  if (!all(is.na(loadingsNOROT))) {				       
    communalities <- as.matrix(diag(loadingsNOROT %*% t(loadingsNOROT))) 	
    communalities <- cbind(smcINITIAL, communalities) 
    rownames(communalities) <- cnoms
    colnames(communalities) <- c('Initial', 'Extraction')
  }
  maxlikeOutput <- list(loadingsNOROT=loadingsNOROT, communalities = communalities)
  
  return(invisible(maxlikeOutput))
}






PA_FA <- function (cormat, Nfactors, Ncases, iterpaf=100) {
  
  # CFA / PAF  (Bernstein p 189; smc = from Bernstein p 104)
  
  Nvars <- dim(cormat)[1]
  
  eigenvalues <- eigen(cormat)$values
  
  cnoms <- colnames(cormat)
  
  converge  <- .001
  rpaf <- as.matrix(cormat)
  smc <- 1 - (1 / diag(solve(rpaf)))
  smcINITIAL <- smc  # initial communalities
  
  for (iter in 1:(iterpaf + 1)) {
    diag(rpaf) <- smc # putting smcs on the main diagonal of r
    eigval <-  diag((eigen(rpaf) $values))
    # substituting zero for negative eigenvalues
    for (luper in 1:nrow(eigval)) { if (eigval[luper,luper] < 0) { eigval[luper,luper] <- 0 }}
    eigvect <- eigen(rpaf) $vectors
    if (Nfactors == 1) {
      loadingsNOROT <- eigvect[,1:Nfactors] * sqrt(eigval[1:Nfactors,1:Nfactors])
      communalities <- loadingsNOROT^2
    }else {
      loadingsNOROT <- eigvect[,1:Nfactors] %*% sqrt(eigval[1:Nfactors,1:Nfactors])
      communalities <- rowSums(loadingsNOROT^2) 
    }
    if (max(max(abs(communalities-smc))) < converge) { break }
    if (max(max(abs(communalities-smc))) >= converge  & iter < iterpaf) { smc <- communalities }
  }
  loadingsNOROT <- as.matrix(loadingsNOROT)
  
  communalities <- as.matrix(diag(loadingsNOROT %*% t(loadingsNOROT))) 	
  communalities <- cbind(smcINITIAL, communalities) 
  rownames(communalities) <- cnoms
  colnames(communalities) <- c('Initial', 'Extraction')
  
  pafOutput <- list(loadingsNOROT=loadingsNOROT, communalities=communalities) 
  
  return(invisible(pafOutput))
}





# Harris, C. W. On factors and factor scores. P 32, 363-379.

# Harris, C. W. (1962). Some Rao-Guttman relationships." Psychometrika, 27,  247-63. 

# HARRIS, CHESTER W. "Canonical Factor Models for the
# Description of Change." Problems in Measuring Change. (Edited by Chester W.
# Harris.) Madison: University of Wisconsin Press, 1963. Chapter 8, pp.
# 138-55. (a) 

# HARRIS, CHESTER W., editor. Problems in Measuring Change. Madison: University of Wisconsin Press, 1963. 259 pp. (b) 

# HARRIS, CHESTER W. "Some Recent Developments in Factor Analysis." Educational and Psychological Measurement 2 4 : 193-206; Summer 1964. 

# HARRIS, CHESTER W., and KAISER, HENRY F. "Oblique Factor Analytic Solutions by Orthogonal Transformations." Psychometrika 29: 347-62; December 1964.

# Guttman, L. (1953). Image theory for the structure of quantitative
# variates. Psychometrika 18, 277-296.




