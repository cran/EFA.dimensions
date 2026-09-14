

# var_partit_mat

# GPA options, LV_options = make sure all elements are provided

# named group_keys matches names in rawdata / cormat?



BIFACTOR <- function(loadings = NULL,
                     rawdata = NULL,   
                     cormat = NULL, Ncases = NULL, corkind = 'pearson', 
                     Nfactors = 4,
                     bifactor_kind = 'bifactorT', 
                     EFA_options    = list(extraction = 'minres', 
                                           rotation = 'oblimin'),
                     LV_options = list(group_keys = NULL, 
                                       estimator='MLR', 
                                       rotation = 'bigeomin',
                                       resid_correls=NULL, 
                                       LV_names=NULL,
                                       ordered = FALSE),
                     schmid_options = list(extraction = 'minres', 
                                           rotation = 'oblimin'),
                     GPA_options = list(delta = .01, 
                                        epsilon = .00001, 
                                        normalize = FALSE, 
                                        maxit = 1000, 
                                        randomStarts = 50),
                     min_loading = .2,
                     verbose = TRUE) { 
  
  
  #############################  argument checks   ##################################
  
  # is the corkind method valid?
  corkind <- corkind_check(corkind)
    
  # is the bifactor_kind valid?
  if (!bifactor_kind %in% c('bifactorQ', 'bifactorT', 'bigeominT', 'bigeominQ', 
                            'SL', 'SLiD', 'DSL', 'CFA', 'ESEM', 'none')) {
    cat('\nThe entry for bifactor_kind, ', bifactor_kind, 
        ', is not one of the options for this function.', sep='')
    cat('\n"bifactorT" will be used instead.')
    bifactor_kind <- 'bifactorT'
  }
  
  # check EFA_options if loadings = NULL
  if (is.null(loadings))  EFA_options <- EFA_options_check(EFA_options)
  
  # if GPArotation will be used, EFA_options$rotation must be 'none'
  if (bifactor_kind %in% c('bifactorQ', 'bifactorT', 'bigeominT', 'bigeominQ')) {
    if (EFA_options$rotation != 'none') {
      cat('\n\nFor bifactor_kind = ', bifactor_kind, 
          ', the EFA rotation must be, and will be changed to, "none"\n', sep='') 
      EFA_options$rotation <- 'none'
    }
  }

  # check schmid_options if schmid will be used
  if (bifactor_kind %in% c('SL', 'SLiD', 'DSL')) 
    schmid_options <- schmid_options_check(schmid_options)
  
  # if schmid will be used, is the EFA_options rotation method valid?
  if (bifactor_kind %in% c('SL', 'SLiD', 'DSL')) {
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
  if (bifactor_kind %in% c('CFA','ESEM'))  LV_options <- LV_options_check(LV_options)
  
  # LV_options$group_keys is required when bifactor_kind = CFA
  if (any(bifactor_kind == 'CFA' & is.null(LV_options$group_keys))) {
    cat('\n\nCFA was specified in bifactor_kind but group_keys in LV_options = NULL')
    cat('\nEither provide LV_options$group_keys or do not request CFA.')
    cat('\nbifactor_kind was changed to bifactorT, to prevent errors')
    bifactor_kind <- 'bifactorT'
  }
    
  # check GPA_options if GPA functions will be used
  if (bifactor_kind %in% c('bifactorQ', 'bifactorT', 'bigeominT', 'bigeominQ'))
    GPA_options <- GPA_options_check(GPA_options)
  

  # rawdata is required for bifactor_kind = CFA & ESEM
  if (is.null(rawdata) & (bifactor_kind == 'CFA' | bifactor_kind == 'ESEM')) 
    stop('\nrawdata is required when bifactor_kind = ', bifactor_kind)
  
  
  
  ############################  rawdata, cormat setup   ##########################
  
  # if both rawdata & cormat are provided, use only rawdata
  if (!is.null(rawdata) & !is.null(cormat)) {
    message('\n\nrawdata & cormat were both provided. cormat will be ignored.')
    cormat <- NULL
  } 
  
  # if bifactor_kind = CFA or ESEM, then the length of group_keys (if provided) must be = # vars in rawdata / cormat
  if (bifactor_kind == 'CFA' | bifactor_kind == 'ESEM') {
    if (!is.null(LV_options$group_keys))  {
      if (!is.null(rawdata) & length(LV_options$group_keys) != ncol(rawdata)) {
        message('\n\nThe number of items in group_keys is not equal to the number of variables in rawdata.')
        message('\nOnly the first ', length(LV_options$group_keys), ' variables in rawdata will be used.')
        rawdata <- rawdata[,1:length(LV_options$group_keys)]
      }
      if (!is.null(cormat)) { 
        if (length(LV_options$group_keys) != ncol(cormat)) {
          message('\n\nThe number of items in group_keys is not equal to the number of variables in cormat.')
          message('\nOnly the first ', length(LV_options$group_keys), ' variables in cormat will be used.')
          cormat <- cormat[,1:length(LV_options$group_keys)]
          cormat <- cormat[1:length(LV_options$group_keys),]
        }
      }
    }
  }
  
  # make sure all vars in rawdata are numbers
  if (!all(sapply(rawdata, is.numeric))) 
    stop('\nNot all of the variables in rawdata are numeric. The analyses cannot proceed.')
 
  if (!is.null(rawdata))  Ncases <- nrow(rawdata)
  
  # cormat, if not provided; required for c('SL', 'SLiD', 'DSL') only, & need rawdata
  if (bifactor_kind %in% c('SL', 'SLiD', 'DSL')) {
    
    if (is.null(cormat) & !is.null(rawdata)) {
      
      cormat <- cor(rawdata)
    }
  }
  
  
  ##################################  Nfactors   ####################################
  
  # if LV_options$group_keys is provided, ignore Nfactors & compute it from LV_options$group_keys
  if (!bifactor_kind %in% c('CFA', 'ESEM') & !is.null(LV_options$group_keys)) 
    Nfactors <- length(unique(LV_options$group_keys))
  
 
  ##############################  loadings_for_BIF   ################################
  
  loadings_for_BIF <- NULL
  
  # no EFA because loadings were provided
  if (!is.null(loadings))  loadings_for_BIF <- loadings 
  
  # run EFA when loadings are needed (for GPArotations only) but not provided 
  # need raw data or cormat & Ncases
  if (bifactor_kind %in% c('bifactorQ', 'bifactorT', 'bigeominT', 'bigeominQ')) {
    
    if (is.null(loadings) & ( (!is.null(cormat) & !is.null(Ncases)) | !is.null(rawdata) )) {
      
       if (!is.null(cormat) & !is.null(Ncases)) {
        efa_output <- EFA(data = cormat, 
                          Nfactors   = Nfactors, 
                          extraction = EFA_options$extraction, 
                          rotation   = EFA_options$rotation, 
                          corkind = corkind, Ncases = Ncases, 
                          iterpaf=100, ppower = 3, 
                          GPA_options = GPA_options, 
                          verbose=FALSE)
      }
      
      if (!is.null(rawdata)) {
        efa_output <- EFA(data = rawdata, 
                          Nfactors   = Nfactors, 
                          extraction = EFA_options$extraction, 
                          rotation   = EFA_options$rotation, 
                          corkind = corkind, 
                          iterpaf=100, ppower = 3, 
                          GPA_options = GPA_options, 
                          verbose=FALSE)
        
        Ncases <- nrow(rawdata)
      }
      
      if (EFA_options$rotation == 'none')                                 
        loadings_for_BIF <- efa_output$loadingsNOROT
      
      if (EFA_options$rotation != 'none' & is.null(efa_output$pattern))   
        loadings_for_BIF <- efa_output$loadingsROT
      
      if (EFA_options$rotation != 'none' & !is.null(efa_output$pattern))  
        loadings_for_BIF <- efa_output$pattern
    }
  }
  
  
  ###########################  bifactor analyses  ###################################
  
  bif_output <- bifactor_engine(loadings = loadings_for_BIF, 
                                cormat = cormat, corkind = corkind, Ncases = Ncases, 
                                rawdata = rawdata,
                                Nfactors = Nfactors,
                                bifactor_kind = bifactor_kind, 
                                schmid_options = schmid_options,
                                LV_options = LV_options,
                                GPA_options = GPA_options,
                                min_loading = min_loading) 
  
  
  #############################  output   ###########################################

    if (verbose)  {
    
    cat('\n\nBIFACTOR output:\n')

    cat('\nThe number of factors = ', Nfactors)
    cat('\nThe number of group factors = ', (Nfactors - 1))
    
    if (is.null(loadings) & bifactor_kind %in% c('bifactorQ', 'bifactorT', 'bigeominT', 'bigeominQ')) {
      cat('\nEFA extraction method = ', EFA_options$extraction)
      cat('\nEFA rotation method   = ', EFA_options$rotation)
    }
    
    if (is.null(loadings) & bifactor_kind %in% c('SL', 'SLiD', 'DSL')) {
      cat('\nschmid extraction method = ', schmid_options$extraction)
      cat('\nschmid rotation method   = ', schmid_options$rotation)
    }
    
    if (!is.null(Ncases)) { cat('\nThe number of cases = ', Ncases, '\n')
    } else { cat('\nThe number of cases is not available\n') }
    
    if (!bifactor_kind %in% c('CFA','ESEM')) {
      cat('\n\nThe loadings that were entered into the bifactor rotations:\n\n')
      if (!is.null(loadings_for_BIF))  print(round(loadings_for_BIF,2))
      if (bifactor_kind %in% c('SL','SLiD','DSL')) print(round(bif_output$loadings_for_BIF,2))
    }
    
    show_bifactor_stats(bif_output)
  }
  
  return(invisible(bif_output))
}

