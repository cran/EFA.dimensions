

# sample_cov = cormat,
# sample.nobs = Ncases,



ESEM <- function(data, 
                 method = 'startvalues', # efa_blocks
                 Nfactors = NULL, 
                 extraction = 'ml', 
                 rotation_EFA = 'oblimin',  
                 rotation_LV  = 'geomin',  
                 corkind = 'pearson', 
                 anchors = NULL, 
                 estimator = 'ML', 
                 ordered = FALSE, 
                 target = NULL, target_keys = NULL, 
                 bifactor = FALSE, 
                 verbose = TRUE) {
  
  if (verbose)  cat('\n\nExploratory Structural Equation Modeling:')
  
  ###############################  data  setup   ####################################
  
  data <- MISSING_DROP(data)
  
  N_obsvd_vars  <- ncol(data)
  
  varnames <- colnames(data)
  
  loadmat <- targ_mat <- NULL
  

  #############################  argument checks   ##################################
  
  # is the corkind method valid?
  corkind <- corkind_check(corkind)
  
  #  extraction & rotation_EFA are checked in the "model syntax for startvalues" section
  
  # check rotation_LV
  if (!rotation_LV %in% c('varimax','quartimax','orthomax','oblimin',
     'quartimin','geomin','promax','entropy','mccammon','infomax',
     'tandem2','oblimax','bentler','simplimax','target','cf','crawford-ferguson',
      'cf-quartimax','cf-varimax','cf-equamax','cf-parsimax','cf-facparsim')) {
    cat('\nThe entry for rotation_LV, ', rotation_LV, 
        ', is not one of the options for this function.', sep='')
    cat('\n"geomin" will be used instead.')
    rotation_LV <- 'geomin'
  }

  # check if Nfactors can be determined
  if (is.null(Nfactors) & is.null(target_keys))  {
    cat('\n\nThe Nfactors & target_keys arguments were both NULL.\n')
    cat('\nThe number of factors for the analyses could not be determined.\n')
    cat('\nNfactors was set to 3 to enable to analyses to proceed.\n')
    Nfactors <- 3
  }

  # check on ordered
  # ordered can be TRUE, FALSE, or it can have the names of ordered variables
  # if ordered != FALSE, check the number of values for each variable & 
  # change to ordered = FALSE if all have > 12 levels (as in the lavaan warning in such cases)
  # change to ordered = the names of the variables with < 13 values, if only some have < 13 values
  ordered <- ordered_data_check(ordered = ordered, rawdata = data[,varnames])
  
  # check if ordered is compatible with estimator & change to WLSMV if not
  estimator <- ordered_estimator_check(ordered = ordered, estimator = estimator)
  
  # # make sure ordered is compatible with estimator
  # if (ordered != FALSE & 
  #     !estimator %in% c('WLSMV','WLSM','DWLS','ULSMV','ULSM','ULS','PML')) {
  #   cat('\n\nThe estimator, ', estimator, ', cannot be used for ordered data.', sep='')
  #   cat('\nThe estimator has therefore been changed to WLSMV.\n')
  #   estimator <- 'WLSMV'
  # }
  
  # # does the estimator work for ordered data?
  # if (ordered & !estimator %in% c('WLSMV', 'DWLS', 'WLS', 'ULS', 'ULSMV')) {
  #   cat('\nThe', estimator, 'estimator does not work for ordered data. WLSMV will be used instead.')
  #   estimator <- 'WLSMV'
  # }
  
  
  #############################  model syntax for efa_blocks   ######################
    
  if (method == 'efa_blocks') {
    
    esem_model_syntax <- lavaan_model(varnames = varnames, 
                                      model = NULL, 
                                      keys = target_keys,
                                      LV_names = NULL, 
                                      resid_correls = NULL,
                                      Nfactors = Nfactors,
                                      bifactor = FALSE,
                                      esem = TRUE)
  }
  
  
  #############################  model syntax for startvalues   #####################
  
  if (method == 'startvalues') {
    
    # is the extraction method valid?
    if (!extraction %in% c('paf', 'ml', 'image', 'minres', 'uls', 'ols', 'wls', 
                           'gls', 'alpha', 'fullinfo')) {
      cat('\n\nThe entry for extraction, ', extraction, ', is not one of the options for this function.', sep='')
      cat('\n"paf" will be used instead.')
      extraction <- 'paf'
    }
    
    # is the rotation_EFA method valid?
    if (!rotation_EFA %in% c('bentlerQ', 'bentlerT', 'bifactorQ', 'bifactorT', 'bigeominQ', 
                         'bigeominT', 'entropy', 'equamax', 'geominQ', 'geominT', 
                         'oblimax', 'oblimin', 'promax', 'quartimax', 'quartimin', 
                         'simplimax', 'varimax', 'targetQ', 'targetT', 'none')) {
      cat('\n\nThe entry for rotation_EFA, ', rotation_EFA, ', is not one of the options for this function.', sep='')
      cat('\n"promax" will be used instead.')
      rotation_EFA <- 'promax'
    }
    
    if (is.null(Nfactors)) {		
      Nfactors <- EMPKC(data=data, corkind=corkind, verbose=FALSE)$NfactorsEMPKC
      NfactorsWasNull <- TRUE
    } else {NfactorsWasNull <- FALSE}
    
    if (!is.null(target_keys)) {
      
      if (!is.null(target)) 
        cat('\n\ntarget_keys were specified but so was target. The target_keys will be ignored')
      
      if (is.null(target)) {
        
        # check for errors in target_keys
        if (length(target_keys) != N_obsvd_vars)
          cat('\n\nThe length of target_keys is not equal to the number of variables in data.')
        if (min(target_keys) != 1)  cat('\nThe smallest number in target_keys is not 1')
        if (!all(seq(min(target_keys):max(target_keys)) %in% unique((target_keys))))
          cat('\n\nNot all of the possible factor numbers appear in target_keys')
        
        target <- matrix(0, length(target_keys), max(target_keys))
        for (lupe in 1:length(target_keys))  target[lupe, target_keys[lupe]] <- NA
        rownames(target) <- varnames
        colnames(target) <- paste('Factor_', 1:ncol(target), sep = '')
      }
    }
    
    if (is.null(target)) {
      
      if (rotation_EFA %in% c('targetQ', 'targetT')) {
        cat('\n\nrotation_EFA was set to', rotation_EFA, 'but a target was provided.')
        cat('\nrotation_EFA was therefore set to "promax" for the EFA.')
        rotation_EFA <- "promax"
      }
      
      efa_output <- EFA(data=data, extraction = extraction, corkind=corkind,  
                        Nfactors=Nfactors, rotation=rotation_EFA, verbose=FALSE)
      
      if (all(!is.na(efa_output$pattern)))  loadmat <- efa_output$pattern
      if (all( is.na(efa_output$pattern)))  loadmat <- efa_output$loadingsROT
    }
    
    if (!is.null(target)) {
      
      if (!rotation_EFA %in% c('targetQ', 'targetT')) {
        cat('\n\nA target was provided and so rotation_EFA was set to "targetQ" for the EFA.')
        rotation_EFA <- "targetQ"
      }
      
      # make sure the size of target is correct
      if (ncol(target) != Nfactors) {
        cat('\nThe number of columns in target', ncol(target), 'was not = Nfactors')
        cat('\nNfactors was therefore changed to the number of columns in target')
        Nfactors <- ncol(target)
      }
      if (nrow(target) != N_obsvd_vars)
        stop('The number of rows in target', ncol(target),
             'was not = to the number of variables in data,', N_obsvd_vars)
      
      # target rotation using psych::fa seems to require that target is a list, not a matrix
      # efa_output <- psych::fa(r = data, nfactors = Nfactors, fm = extraction,
      #                         rotate = rotation_EFA, Target = target)
      #                         # scores=scores,
      #                         # residuals=residuals,
      #                         # missing=missing)
      # 
      # loadmat <- efa_output$loadings
      
      # first run an efa with no rotation
      efa_output <- EFA(data=data, extraction = extraction, corkind=corkind,  
                        Nfactors=Nfactors, rotation='none', verbose=FALSE)$loadingsNOROT
      
      # use GPArotation for target rotation
      if (rotation_EFA == "targetQ")  
        loadmat <- targetQ(efa_output, Target = target)$loadings[1:nrow(efa_output),]
      if (rotation_EFA == "targetT")  
        loadmat <- targetT(efa_output, Target = target)$loadings[1:nrow(efa_output),]
    }
    
    
    colnames(loadmat) <- paste('Factor_', 1:ncol(loadmat), sep = '')
    
    
    # # loadmat <- zapsmall(matrix(round(efa_output, 2), nrow = N_obsvd_vars, ncol = Nfactors))
    # loadmat <- zapsmall(matrix(efa_output, nrow = N_obsvd_vars, ncol = Nfactors))
    # rownames(loadmat) <- varnames
    
    
    # if anchors are provided, test if they are valid rownames
    if (!is.null(anchors)) {
      if (!all(anchors %in% varnames)) 
        stop('anchors were specified but they do not all appear in the variable names')
    }
    # if anchors are not provided, id them - the highest-loading item on each factor
    if (is.null(anchors)) 
      anchors <- rownames(loadmat)[apply(loadmat, 2, function(x) which(abs(x) == max(abs(x))))]
    
    # create a set of lavaan equations from the efa loadings
    # use the efa loadings as start values, except do not use "start" for the anchor items
    terms  <- vector()
    for (lupe in 1:Nfactors) {
      
      dum <- paste0("F", lupe,"=~")
      
      for (lupe_rows in 1:nrow(loadmat)) {
        
        if (lupe_rows != nrow(loadmat)) {
          
          if (rownames(loadmat)[lupe_rows] %in% anchors)  dum <- 
              paste0(dum, paste0(loadmat[lupe_rows,lupe], "*", rownames(loadmat)[lupe_rows], "+\n"))
          if (!rownames(loadmat)[lupe_rows] %in% anchors)  dum <- 
              paste0(dum, paste0('start(', loadmat[lupe_rows,lupe], ")*", rownames(loadmat)[lupe_rows], "+\n"))
        } 
        
        if (lupe_rows == nrow(loadmat)) {
          
          if (rownames(loadmat)[lupe_rows] %in% anchors)  dum <- 
              paste0(dum, paste0(loadmat[lupe_rows,lupe], "*", rownames(loadmat)[lupe_rows], "\n"))
          if (!rownames(loadmat)[lupe_rows] %in% anchors)  dum <- 
              paste0(dum, paste0('start(', loadmat[lupe_rows,lupe], ")*", rownames(loadmat)[lupe_rows], "\n"))
        } 
      }     
      terms[lupe] <- dum
    }
    
    if (bifactor) {
      
      dum_gen <- paste0("Gen=~")   
      
      for (lupe_rows in 1:nrow(loadmat)) {
        
        if (lupe_rows != nrow(loadmat)) 
          
          dum_gen <- paste0(dum_gen, rownames(loadmat)[lupe_rows], "+\n")
        
        if (lupe_rows == nrow(loadmat)) 
          
          dum_gen <- paste0(dum_gen, rownames(loadmat)[lupe_rows], "\n")
      }
      # dum_gen
      # writeLines(dum_gen)
      terms <- c(dum_gen, terms)
      # terms
      # writeLines(terms)
    }
    
    esem_model_syntax <- paste(terms, collapse = "\n")
  }
  
  
  ###################################  lavaan   #####################################
  
  # cat('\n\nesem_model_syntax\n'); writeLines(esem_model_syntax)
  
  # Test the syntax validity
  tryCatch({
    parsed_table <- lavaanify(model = esem_model_syntax)
    # message('The lavaan model syntax is valid.')
    # message('\nLook at the parameter table structure:\n')
    # print((parsed_table[, c("lhs", "op", "rhs")]))
  }, error = function(e) {
    message("Syntax Error found:")
    print(e$message)
  })
  
  
  # # "Correlated uniqueness"
  # terms[6] <- "A1 ~~ C2+E3+N3\n C2 ~~ E3+N3\n E3 ~~ N3"
  
  
  if (method == 'efa_blocks' & !is.null(target_keys)) {
    
    # generate a target matrix from keys
    targ_mat <- target_matrix(keys=target_keys, varnames, bifactor=FALSE)
    
    rotation_LV <- "target"
    
    lavaan_output <- lavaan::cfa(esem_model_syntax, data=data, 
                                 estimator=estimator, #  "WLSMV", 
                                 std.lv=TRUE, 
                                 ordered=ordered, 
                                 rotation = rotation_LV, 
                                 rotation.args = list(target = targ_mat),
                                 verbose=FALSE)
  }
  
  
  if (method == 'startvalues' | (method == 'efa_blocks' & is.null(target_keys))) {
    
    lavaan_output <- lavaan::cfa(esem_model_syntax, data=data, 
                                 estimator=estimator, #  "WLSMV", 
                                 rotation = rotation_LV,
                                 std.lv=TRUE, 
                                 ordered=ordered, 
                                 verbose=FALSE)
  }
  
  # fits <- lavaan::fitmeasures(lavaan_output, c("cfi.robust","tli.robust","rmsea.robust","srmr"))
  # print(fits)
  
  
  #############################  output   ###########################################

    if (verbose) {
    
    cat('\n\nThe number of useable cases in data:', nrow(data))
    
    cat('\n\nThe number of variables:', N_obsvd_vars)
    
    cat('\n\nmethod = ', method)
    
    if (method == 'startvalues') {
      
      cat('\n\nThe efa extraction method:', extraction)
      
      if (!rotation_EFA %in% c('targetQ','targetT')) 
        cat('\n\nThe efa rotation_EFA method:', rotation_EFA)
      
      if (rotation_EFA == 'targetQ')  
        cat('\n\nThe efa rotation_EFA method:', rotation_EFA, ' (oblique)')
      
      if (rotation_EFA == 'targetT')  
        cat('\n\nThe efa rotation_EFA method:', rotation_EFA, ' (orthogonal)')
      
      cat('\n\nThe kind of correlations for the efa:', corkind,'\n')
      
      if (NfactorsWasNull ) {
        cat('\nNfactors was not specified and so the EMPKC test was')
        cat('\nconducted to determine the number of factors to extract: Nfactors = ', Nfactors)		
      } else if (!NfactorsWasNull) {
        cat('\nThe number of factors: ', Nfactors)
      }
      
      if (!is.null(target)) {
        cat('\n\nThe target loading matrix:\n\n')
        print(round(target,3))
      }
      
      cat('\n\nThe efa loadings:\n\n')
      print(round(loadmat,3))
      
      cat('\n\nThe identified anchors: ', anchors, '\n')
    }
    
    if (!is.null(targ_mat)) {
      cat('\n\nThe target loading matrix:\n\n')
      print(round(targ_mat,3))
    }
    
    # cat('\n\nESEM results: \n')
    
    cat('\n\nThe lavaan rotation method (rotation_LV):', rotation_LV, '\n')
    
    show_lavaan_stats(lavaan_output, these = 'all')
    
    # print(lavaan::summary(lavaan_output, fit.measures=TRUE, standardized=TRUE))
    # lavaanPlot(model = lavaan_output, coefs = TRUE,
    #            stand = TRUE,
    #            edge_options = list(color ='grey'))
    
  }
  
  output <- list(loadings = loadmat, esem_model_syntax = esem_model_syntax, 
                 anchors = anchors, targ_mat = targ_mat,
                 lavaan_output = lavaan_output)
  
  return(invisible(output))
}

