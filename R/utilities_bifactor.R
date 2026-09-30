


bifactor_engine <- function(loadings = NULL, 
                            cormat = NULL, corkind = 'pearson', Ncases = NULL, 
                            rawdata = NULL,
                            Nfactors = 4,
                            bifactor_kind, 
                            schmid_options = list(extraction = 'minres', 
                                                  rotation = 'oblimin'),
                            LV_options = list(group_keys = NULL, 
                                              estimator = 'MLR', 
                                              rotation = 'bigeomin',
                                              resid_correls=NULL, 
                                              LV_names=NULL,
                                              ordered = FALSE),
                            GPA_options = list(delta = .01, 
                                               epsilon = .00001, 
                                               normalize = FALSE, 
                                               maxit = 1000, 
                                               randomStarts = 50), 
                            min_loading = .2) { 
  
  # Nfactors should be the # of group factors + 1
  
  
  #############################  argument checks   ##################################
  
  # just one bifactorkind
  if (length(bifactor_kind) > 1) {
    cat('\nbifactor_kind can have only one value. It has been changed to ', bifactor_kind[1])
    bifactor_kind <- bifactor_kind[1]
  }
  
  
  #############################  set up data & cormat   #############################

  # cormat & Ncases, if not provided & rawdata is provided
  if (is.null(cormat) & !is.null(rawdata)) {
    cormat <- cor(rawdata)
    Ncases = nrow(rawdata)
  }
  
  N_group_factors <- Nfactors - 1
  
  ###########################  bifactor analyses  ###################################
  
  loadings_for_BIF <- loadings
  
  structure_BIF <- phi_BIF <- varexplBIF <- lav_model <- lavaan_output <- NULL
  
  # when loadings are already bifactor loadings
  if (bifactor_kind == 'none')  loadings_BIF <- loadings 
  
  if (bifactor_kind == 'CFA' | bifactor_kind == 'ESEM') {
    
    varnames <- colnames(rawdata)
    
    # check on ordered
    # ordered can be TRUE, FALSE, or it can have the names of ordered variables
    # if ordered != FALSE, check the number of values for each variable & 
    # change to ordered = FALSE if all have > 10 levels 
    # change to ordered = the names of the variables with <= 10 values, if only some have > 10 values
    LV_options$ordered <- ordered_data_check(ordered = LV_options$ordered, rawdata = rawdata)
    
    # check if ordered is compatible with estimator & change to WLSMV if not
    LV_options$estimator <- ordered_estimator_check(ordered = LV_options$ordered, 
                                                    estimator = LV_options$estimator)
    
    if (bifactor_kind == 'CFA') {
      
      # generate lavaan CFA model syntax
      lav_model <- lavaan_model(varnames = varnames, 
                                model=NULL, 
                                keys = LV_options$group_keys,  
                                Nfactors = N_group_factors,
                                LV_names=NULL, 
                                resid_correls=NULL,
                                bifactor = TRUE,
                                esem = FALSE)
      
      # lavaan CFA
      lavaan_output <- cfa(model = lav_model, 
                           data = rawdata,
                           # sample_cov = cormat,
                           # sample.nobs = Ncases,
                           ordered = LV_options$ordered, 
                           estimator = LV_options$estimator, 
                           orthogonal = TRUE)
    }
    
    if (bifactor_kind == 'ESEM') {
      
      # generate lavaan B-ESEM model syntax
      lav_model <- lavaan_model(varnames = colnames(rawdata), 
                                model = NULL, 
                                keys = LV_options$group_keys, 
                                Nfactors = N_group_factors,
                                LV_names = NULL, 
                                resid_correls = NULL,
                                bifactor = TRUE,
                                esem = TRUE)
      
      if (!is.null(LV_options$group_keys)) {
        target <- target_matrix(keys = LV_options$group_keys, varnames=varnames, bifactor=TRUE)
        
        lavaan_output <- lavaan::cfa(lav_model, data=rawdata, 
                                     estimator = LV_options$estimator,
                                     ordered = LV_options$ordered,
                                     std.lv=TRUE,  
                                     rotation = "target", 
                                     rotation.args = list(target     = target, 
                                                          orthogonal = TRUE,
                                                          rstarts    = as.integer(5L)),
                                     verbose=FALSE)
      }
      
      if (is.null(LV_options$group_keys)) {
        lavaan_output <- lavaan::cfa(lav_model, 
                                     data=rawdata, 
                                     ordered = LV_options$ordered,
                                     estimator = LV_options$estimator,
                                     rotation = LV_options$rotation,
                                     # orthogonal = TRUE,
                                     std.lv=TRUE,  
                                     verbose=FALSE)
      }
    }
    
    loadings_BIF <- lavInspect(lavaan_output, what = "std")$lambda
    phi_BIF      <- lavInspect(lavaan_output, what = "std")$psi
    class(loadings_BIF) <- class(phi_BIF) <- "matrix"
  }

  if (bifactor_kind %in% c('bifactorT', 'bifactorQ', 'bigeominT', 'bigeominQ')) {
    
    outp_rotation <- rotation_func(rotation=bifactor_kind, 
                                   loadingsNOROT = loadings_for_BIF, 
                                   GPA_options,
                                   ppower = 3)
    
    if (!is.null(outp_rotation$loadingsROT)) loadings_BIF  <- outp_rotation$loadingsROT
    if (!is.null(outp_rotation$pattern))     loadings_BIF  <- outp_rotation$pattern
    if (!is.null(outp_rotation$structure))   structure_BIF <- outp_rotation$structure
    if (!is.null(outp_rotation$phi))         phi_BIF       <- outp_rotation$phi
  }
  
  if (bifactor_kind == 'SL') {
    
    # Nfactors for bifactor analyses:
    # SL, SLiD, & DSL all do their own factor extractions in which Nfactors = N_group_factors
    # & the general factor is later added to the Nfactors loading matrix in these cases.
    
    schmid_outp <- suppressMessages(
      psych::schmid(model = cormat, 
                    nfactors = N_group_factors, 
                    fm = schmid_options$extraction, 
                    rotate = schmid_options$rotation, 
                    digits=2, n.obs=Ncases, option="equal",
                    Phi=NULL, covar=FALSE, two.ok=FALSE))  #, plot=FALSE)
    
    loadings_SL <- schmid_outp$sl
    
    loadings_BIF <- loadings_SL[, -which(colnames(loadings_SL) %in% c('h2','u2','p2','com'))]
    
    loadings_for_BIF <- schmid_outp$oblique
  }
  
  if (bifactor_kind == 'SLiD') {
    
    # SLiD_rotation - 2021 Garcia-Garzon - On Omega Hierarchical Estimation - 
    # A Comparison of Exploratory Bi-Factor Analysis Algorithms
    
    SLiD_output <- suppressMessages(
      slid_rotation(data = cormat, 
                    n_factors = N_group_factors, 
                    fm = schmid_options$extraction, 
                    rotate = schmid_options$rotation, 
                    max_iter = 100, tol = 1e-5))
    
    loadings_BIF <- SLiD_output$bifactor_loadings
    
    loadings_for_BIF <- SLiD_output$loadings_for_BIF
  }
  
  if (bifactor_kind == 'DSL') {  
    
    # Direct Schmid Leiman -- psych
    directSl_outp <- suppressMessages(
      psych::directSl(cormat, 
                      nfactors = N_group_factors, 
                      fm = schmid_options$extraction, 
                      rotate = schmid_options$rotation, cut=.3))  #$direct  # $f$loadings
    
    loadings_DSL <- unclass(directSl_outp$direct)
    
    loadings_BIF <- loadings_DSL
    
    loadings_for_BIF <- directSl_outp$f$loadings[,1:N_group_factors]
  }
  
  
  if (anyNA(loadings_BIF))
    cat('\n\nSome loadings could not be estimated and are NA. Expect problems.\n\n')
  
  colnames(loadings_BIF) <- 
    c('General', c(paste('Group', 1:(ncol(loadings_BIF) - 1), sep=' ')))
  
  if (is.null(rownames(loadings_BIF))) 
    rownames(loadings_BIF) <- paste('v', 1:nrow(loadings_BIF), sep='')
  
  if (!is.null(structure_BIF))  {
    
    colnames(structure_BIF) <- 
      c('General', c(paste('Group', 1:(ncol(structure_BIF) - 1), sep=' ')))
    
    if (is.null(rownames(structure_BIF))) 
      rownames(structure_BIF) <- paste('v', 1:nrow(structure_BIF), sep='')
  }
  
  if (!is.null(phi_BIF))  {
    colnames(phi_BIF) <- rownames(phi_BIF) <- 
    c('General', c(paste('Group', 1:(ncol(phi_BIF) - 1), sep=' ')))
  }
  
  # check if loadings are consistent with a bifactor structure, based on min_loadings
  # i.e., if any item has > 2 loadings that are >= min_loading 
  # re: an item should load on the general factor & on just one group factor
  
  bifactor_flag <- FALSE
  if (any(rowSums(abs(loadings_BIF) >= min_loading) > 2)) { bifactor_flag <- TRUE }    
  
  
  ###########################  omega & bifactor stats  ##############################
  
  omega_total <- omega_hierl <- Calpha.z <-ECV <- ARPB <- FD <- PUC <- coef_H <- 
    rmsr <- rmsr_gen <- rmsr_psych <- var_partit_mat <- NULL
  
  # omega
  omega_total <- omega_t(loadings_BIF, min_loading = min_loading) 
  omega_hierl <- omega_h(loadings_BIF, min_loading = min_loading)
  # omega_total <- omega_t(loadings_BIF, cormat) 
  # omega_hierl <- omega_h(loadings_BIF, cormat)
  
  # standardized alpha
  if (!is.null(cormat))  Calpha.z <- Cronbach.alpha.z(cormat)
  
  
  # ECV - Explained Common Variance of the general & specific factors (3 versions, SS, SG, GS) 

  # # ECV - Explained Common Variance of the (first) general factor
  # eigenvalues <- colSums(loadings_BIF^2) 
  # ECV <- eigenvalues[1] / sum(eigenvalues)

  # putting the general factor ECV as the first value in all 3
  # 2023 Dueber, Toland - A Bifactor Approach to Subscore Assessment  p 224 formula 3
  # ECV_SS <- ECV_SG <- ECV_GS <- 
  #   (sum(loadings_BIF[,1]^2) / (sum(loadings_BIF[,1]^2) + 
  #                                 sum(colSums(loadings_BIF[,2:ncol(loadings_BIF), drop=FALSE]^2))))
  
  ECV_SS <- ECV_SG <- ECV_GS <- NULL
  
  for (lupe in 1:ncol(loadings_BIF)) {
    
    # extract loadings based on min_loading
    group_TF <- abs(loadings_BIF[,lupe]) >= min_loading
    # dum <- cbind( loadings_BIF[group_TF,1], loadings_BIF[group_TF,lupe])
    # dum <- loadings_BIF[group_TF,]
    dum <- loadings_BIF[group_TF, , drop = FALSE]   # thanks Yves
    
    # ECV_SG (Specific-dimension Explained Common Variance, also called ECV S&E) is 
    # the proportion of the total common variance across all items that is explained 
    # by a specific (group) factor rather than the general factor.     # ECV_SG <- c(ECV_SG, (sum(dum[,2]^2) / ( sum(loadings_BIF^2))))
    ECV_SG <- c(ECV_SG, (sum(dum[,lupe]^2) / ( sum(loadings_BIF^2))))
    
    # ECV_SS (Explained Common Variance of a Specific Factor with respect to itself) 
    # is the proportion of common variance for the items in a specific subscale/group 
    # factor that is explained by that specific factor relative to the total common 
    # variance (from both general and specific factors) captured by those same items. 
    # 2023 Dueber, Toland - A Bifactor Approach to Subscore Assessment  p 224 formula 6
    # ECV_SS <- c(ECV_SS, ( sum(dum[,2]^2) / ( sum(dum[,1]^2) + sum(dum[,2]^2)) ) )
    ECV_SS <- c(ECV_SS, ( sum(dum[,lupe]^2) / ( sum(dum^2) ) ) )
    
    # ECV_GS (General-Specific Explained Common Variance, or the within-domain ECV 
    # for a specific factor) is the proportion of common variance in the 
    # indicators of a specific subscale/group factor that is driven by the general 
    # factor.
    # ECV_GS <- c(ECV_GS, (sum(dum[,1]^2) / ( sum(dum[,1]^2) + sum(dum[,2]^2))))
    # ECV_GS <- c(ECV_GS, (sum(dum[,1]^2) / ( sum(dum[,1]^2) + sum(dum[,lupe]^2))))
    ECV_GS <- c(ECV_GS, (sum(dum[,1]^2) / ( sum(dum^2) )))
  }

  # bifactorIndices also computes the above 3 ECV coefficients, but not for ECV-GS
  # or ECV_SG when there are no zero loadings in a column (i.e., when all items
  # have non-zero loadings on all factors)
  
    
  # PUC    not using, it is based on counts of 0-value loadings, as in CFA
  # the number of items on each factor that are >= min_loading
  # if any item has > 2 loadings that are >= min_loading, then do not compute PUC 
  if (bifactor_flag) { PUC <- NA
  } else {
    N_items_ge_min <- colSums( abs(loadings_BIF) >= min_loading)
    Nitems <- nrow(loadings_BIF)
    N_correls_tot <- (Nitems*(Nitems-1) / 2)
    N_correls_contaminated <- sum(N_items_ge_min * (N_items_ge_min - 1) / 2) -
      N_correls_tot
    PUC <- 1 - N_correls_contaminated / N_correls_tot
  }
  
  # ARPB -- Average Relative Parameter Bias (& and Relative Parameter Bias, for items)
  if (!is.null(cormat)) {
    
    # needs unidimensional model loadings 
    if (!bifactor_kind %in% c('CFA','ESEM'))
      loadings_BIF_unid <- EFA.dimensions::EFA(data=cormat, 
                                               extraction = schmid_options$extraction,  
                                               corkind=corkind, Ncases=Ncases, Nfactors = 1, 
                                               rotation='none', verbose=FALSE)$loadingsNOROT
    if (bifactor_kind %in% c('CFA','ESEM')) {
      varnames <- colnames(rawdata)
      # generate lavaan CFA model syntax
      lav_mod_unid  <- vector()
      dum <- paste(paste0('Gen', " =~ "), paste(varnames, collapse = ' + '))
      lav_mod_unid <- paste(lav_mod_unid, '\n', dum, collapse = '\n', sep='')
      CFA_bif_mod_unid <- cfa(model = lav_mod_unid, 
                              data = rawdata,
                              ordered = LV_options$ordered,
                              estimator = LV_options$estimator, # MLR
                              orthogonal = TRUE)
      loadings_BIF_unid <-  lavInspect(CFA_bif_mod_unid, what = "std")$lambda
      class(loadings_BIF_unid) <-  "matrix"
    }
    
    loadings_BIF_gen <- loadings_BIF[,1]
    
    ARPB_items <- abs((loadings_BIF_unid - loadings_BIF_gen) / loadings_BIF_gen)
    
    ARPB_total <- mean(ARPB_items)
    
    ARPB <- list(ARPB_items = ARPB_items, ARPB_total = ARPB_total)
  }
  
  # FD -- factor determinacy index
  # phi is a matrix of factor intercorrelations. For the bifactor model,
  # this matrix always will have ones on the diagonal and zeros elsewhere
  phi_BIF_FD <- diag(1, nrow = ncol(loadings_BIF))
  cormat_reprod <- reproduced_R(loadings_BIF)
  FD <- sqrt(diag(phi_BIF_FD %*% t(loadings_BIF) %*% solve(cormat_reprod) %*% loadings_BIF %*% phi_BIF_FD))
  
  
  # coefficient H -- using the general factor loadings
  coef_H <- coefficient_H(loadings=loadings_BIF)
  # H <-  1/(1+1/(colSums(loadings_BIF^2/(1-loadings_BIF^2))))
  
  
  # rmsr
  if (!is.null(cormat)) {
    
    cormat_reproduced <- reproduced_R(loadings_BIF)
    
    # rmsr for the factor model (general + group factors)
    rmsr <- RMSR_boc(cormat, cormat_reproduced)
    # residuals <- cormat - cormat_reproduced
    # residuals.upper <- as.matrix(residuals[upper.tri(residuals, diag = FALSE)])
    # rmsr <- sqrt(mean(residuals.upper^2)) # rmr is perhaps the more common term for this stat
    
    # rmsr for the general factor-only model
    rmsr_gen <- RMSR_boc(cormat, reproduced_R(loadings_BIF[,1, drop=FALSE]))
    
    # rmsr - psych
    residuals <- cormat - cormat_reproduced
    rstar.off <- sum(residuals^2)/2
    n <- nrow(loadings_BIF)  #number of variables
    rmsr_psych <- sqrt(rstar.off/(n*(n-1)))  #this is the empirical rmsea
  }
  
  
  # item_stats for orthogonal loadings
  if (!bifactor_kind %in% c('bifactorQ', 'bigeominQ')) {
    
    loadings_BIF_sqd <- loadings_BIF^2
    
    h2 <- rowSums(loadings_BIF_sqd)
    
    u2 <- 1 - h2
    
    IECV <- loadings_BIF_sqd[,1] / h2
  }
  
  # item_stats for oblique (pattern) loadings
  if (bifactor_kind %in% c('bifactorQ', 'bigeominQ')) {
    
    # structure_BIF <- loadings_BIF %*% phi_BIF
    
    pattern_structure <- loadings_BIF * structure_BIF
    
    h2 <- rowSums(pattern_structure)
    
    u2 <- 1 - h2
    
    IECV <- pattern_structure[,1]^2 / h2
  }
  
  if (!is.null(cormat)) {
    item_stats <- cbind(h2, u2, IECV, ARPB_items)
    colnames(item_stats) <- c('Communalities', 'Uniquenesses', 'IECV', 'Rel. Param. Bias')
  } else {
    item_stats <- cbind(h2, u2, IECV)
    colnames(item_stats) <- c('Communalities', 'Uniquenesses', 'IECV')
  }
  
  # variance partitioning -- from "GPA3bifactor - Bernaards.pdf"
  loadings_BIF_gen    <- loadings_BIF[, 1]
  loadings_BIF_group  <- loadings_BIF[, -1, drop=FALSE]
  theta <- 1 - rowSums(loadings_BIF^2)
  denom <- sum(loadings_BIF_gen)^2 + sum(loadings_BIF_group^2) + sum(theta)
  
  var_partit_mat <- sum(loadings_BIF_gen)^2 / denom
  for (j in 1:ncol(loadings_BIF_group)) 
    var_partit_mat <- rbind(var_partit_mat, (sum(loadings_BIF_group[, j]^2) / denom) )
  var_partit_mat <- rbind(var_partit_mat, (sum(theta) / denom) )
  
  for (lupe in 1:nrow(var_partit_mat)) {
    
    if (lupe == 1) rownames(var_partit_mat)[1] <- 'General factor   '
    
    if (lupe > 1 & lupe < nrow(var_partit_mat))
      rownames(var_partit_mat)[lupe] <- paste('Group factor', (lupe-1), '  ')
    
    if (lupe == nrow(var_partit_mat)) 
      rownames(var_partit_mat)[lupe] <- 'Measurement error'
  }
  
  if (!is.null(cormat))
    varexplBIF <- VarianceExplained(eigenvalues = eigen(cormat)$values, 
                                    loadingsROT = loadings_BIF, phi = phi_BIF)
  
  
  ###############################  output  ##########################################
  
  output <- list(omega_total = omega_total, omega_hierl = omega_hierl,
                 Calpha.z = Calpha.z,
                 ECV_SS = ECV_SS, ECV_GS = ECV_GS, ECV_SG = ECV_SG,
                 ARPB = ARPB, FD = FD, coef_H = coef_H, PUC = PUC,
                 rmsr = rmsr, rmsr_gen = rmsr_gen, rmsr_psych = rmsr_psych, 
                 var_partit_mat = var_partit_mat, 
                 loadings_BIF = loadings_BIF, loadings_for_BIF = loadings_for_BIF,
                 structure_BIF = structure_BIF, 
                 phi_BIF = phi_BIF, varexplBIF = varexplBIF, 
                 item_stats = item_stats, bifactor_kind = bifactor_kind,
                 min_loading = min_loading, lav_model=lav_model,
                 lavaan_output = lavaan_output, LV_options = LV_options,
                 GPA_options = GPA_options,
                 bifactor_flag = bifactor_flag    
  )
  
  return(invisible(output))
}




show_bifactor_stats <- function(bifactor_output) {
  
  cat('\nbifactor_kind = ', bifactor_output$bifactor_kind)
  
  if (bifactor_output$bifactor_kind %in% 
      c('bifactorT', 'bifactorQ', 'bigeominT', 'bigeominQ')) {
    cat('\n\nGPA_options:')
    cat('\n   delta = ', bifactor_output$GPA_options$delta)
    cat('\n   epsilon = ', bifactor_output$GPA_options$epsilon)
    cat('\n   normalize = ', bifactor_output$GPA_options$normalize)
    cat('\n   maxit = ', bifactor_output$GPA_options$maxit)
    cat('\n   randomStarts = ', bifactor_output$GPA_options$randomStarts)
  }
  
  # cat('\n\nBifactor results:')
  
  if (bifactor_output$bifactor_kind %in% c('CFA','ESEM')) {
    
    if (!is.null(bifactor_output$LV_options$group_keys)) {
      if (is.null(names(bifactor_output$LV_options$group_keys)))
        names(bifactor_output$LV_options$group_keys) <- rownames(bifactor_output$item_stats)
      
      cat('\n\nLV_options group_keys:\n\n'); print(bifactor_output$LV_options$group_keys)
    }
    
    if (is.null(bifactor_output$LV_options$group_keys))
      cat('\n\nLV_options group_keys = none provided') 
      
    cat('\n\nLV_options estimator = ', bifactor_output$LV_options$estimator)
    
    if (is.null(bifactor_output$LV_options$group_keys))
      cat('\n\nLV_options rotation = ', bifactor_output$LV_options$rotation)
    
    cat('\n\nLV_options ordered = ', bifactor_output$LV_options$ordered)
    
    if (!is.null(bifactor_output$LV_options$resid_correls)) {
      cat('\n\nLV_options resid_correls = ', bifactor_output$LV_options$resid_correls)
    } # else { cat('\n\nLV_options resid_correls = none provided') }
    
    cat('\n\nlavaan model syntax:\n\n')
    writeLines(bifactor_output$lav_model)
    
    show_lavaan_stats(bifactor_output$lavaan_output, these = 'fits')
  }
  
  if (bifactor_output$bifactor_flag) {
    cat('\n\n\nWARNING: The bifactor loadings do not seem consistent with a bifactor structure:')
    cat('\n         There are items that have more than two loadings that are >= the')
    cat('\n         min_loading value of', bifactor_output$min_loading)
  }  
  
  if (is.null(bifactor_output$structure_BIF)) {
    cat('\n\n\nBifactor Loadings\n\n')	
    print(round(bifactor_output$loadings_BIF,2), print.gap=3) 
    
    if (bifactor_output$bifactor_kind %in% c('CFA','ESEM')) {
      cat('\nFactor correlations\n\n')
      print(round(bifactor_output$phi_BIF,2), print.gap=4)
    }
  }
  
  if (!is.null(bifactor_output$structure_BIF)) { 
    
    cat('\n\n\n', paste(bifactor_output$bifactor_kind, 
                        'Pattern Matrix (standardized factor loadings)'), '\n\n')	
    print(round(bifactor_output$loadings_BIF,2), print.gap=3)
    
    cat('\n\n\n', paste(bifactor_output$bifactor_kind, 'Structure Matrix'), '\n\n')	
    print(round(bifactor_output$structure_BIF,2), print.gap=3)
    
    cat('\n\n\n', paste(bifactor_output$bifactor_kind, 'Factor Correlations'), '\n\n')	
    print(round(bifactor_output$phi_BIF,2), print.gap=3)
  }
  
  cat('\n\nBifactor Model Statistics')	
  
  cat('\n\n   Omega total =', round(bifactor_output$omega_total[1],2),
      '\n\n   Omega hierarchical =', round(bifactor_output$omega_hierl[1],2) )
  
  if (!is.null(bifactor_output$Calpha.z))
    cat('\n\n   Cronbach alpha =', round(bifactor_output$Calpha.z,2))
  
  if (!is.null(bifactor_output$ARPB$ARPB_total))
    cat('\n\n   Average relative parameter bias (ARPB) =', round(bifactor_output$ARPB$ARPB_total,2))
  
  if (!is.null(bifactor_output$PUC))
    cat('\n\n   Percent of uncontaminated correlations (PUC) =', round(bifactor_output$PUC,2))
  
  if (!is.null(bifactor_output$rmsr))
    cat('\n\n   Root mean square of the residuals using the general & group factors =', round(bifactor_output$rmsr,2))
  
  if (!is.null(bifactor_output$rmsr_gen))
    cat('\n\n   Root mean square of the residuals using only the general factor =', round(bifactor_output$rmsr_gen,2))
  
  
  cat('\n\n\n(min_loading = ', round(bifactor_output$min_loading,4), 
      ', which is important for PUC and for omega and ECV group factor statistics)\n', sep='')
  
  
  cat('\n\nBifactor Factor Statistics\n')	
  factor_stats <- rbind(bifactor_output$omega_total, bifactor_output$omega_hierl, 
                        bifactor_output$ECV_SS, bifactor_output$ECV_SG, bifactor_output$ECV_GS, 
                        bifactor_output$FD, bifactor_output$coef_H)
  rownames(factor_stats) <- c('Omega total', 'Omega hierarchical',
                              'ECV - SS', 'ECV - SG', 'ECV - GS', 
                              'Factor Determinacy', 'coefficient H')
  colnames(factor_stats) <- c('General', c(paste('Group', 1:(ncol(bifactor_output$loadings_BIF) - 1), sep=' ')))
  writeLines(paste0("   ", capture.output(print(round(factor_stats,2), print.gap=3))))
  
  
  cat('\n\nBifactor Item Statistics\n\n')	
  print(round(bifactor_output$item_stats,3), print.gap=4)
  
  
  # cat('\n\nBifactor Partitions of the Total Score Variance\n')	
  # var_partit_mat <- round(bifactor_output$var_partit_mat, 3)
  # for (lupe in 1:nrow(var_partit_mat)) 
  #   cat("\n   ", rownames(var_partit_mat)[lupe], '  ', var_partit_mat[lupe,])
  # 
  # cat('\n\n    The variance partition results indicate that ', var_partit_mat[1,]*100,
  #     '% of the total score', sep = '')
  # cat('\n    variance could be attributable to the general factor, with the group')
  # cat('\n    factors together contributing ',
  #     sum( c(var_partit_mat[2:(nrow(var_partit_mat)-1),]))*100,
  #     '% and measurement error ', var_partit_mat[nrow(var_partit_mat),]*100, '%.', sep = '')
  
  #  The variance explained by each factor may be informative: the general factor accounts for
  # approximately 56.3% of item variance, while the two group factors account for 6.5% and 7.9%
  # respectively, for a total of 70.7%. 
  
  if (!is.null(bifactor_output$varexplBIF)) {
    cat('\n\n\nEigenvalues and Proportions of Total Variance Explained\n')
    cat('\n               Initial            Bifactor\n')  
    print(bifactor_output$varexplBIF, print.gap=2)
  }
  
}  





# coefficient H is a measurement of construct replicability and is
# defined as the extent to which a set of items represents a latent
# variable (Hancock & Mueller, 2001)
coefficient_H <- function(loadings=NULL, cormat=NULL, corkind='pearson',
                          Ncases=NULL, extraction=NULL) {
  
  if (is.null(loadings)) {
    
    efa_output <- EFA.dimensions::EFA(data=cormat, extraction = extraction, 
                                      corkind=corkind, Ncases=Ncases, Nfactors = 1, 
                                      rotation='none', verbose=FALSE)
    
    loadings <- efa_output$loadingsNOROT
    # errors <- 1 - efa_output$communalities
  }
  
  # 2017 McNeish p 418
  proportions <- (loadings^2) / (1 - loadings^2)
  if (is.vector(proportions)) dim(proportions) <- c(length(proportions), 1)
  # sum_of_props <- sum(proportions)
  sum_of_props <- colSums(proportions)
  coef_H <- sum_of_props / (1 + sum_of_props)
  
  # # reliable package approach
  # Lambda <- loadings
  # Psi <- diag(as.vector(errors))
  # numerator <- crossprod(Lambda * 1 / diag(Psi), Lambda)
  # H <- unname( numerator / (1 + numerator) )
  
  # H <-  1/(1+1/(colSums(loadings^2/(1-loadings^2))))
  
  return(invisible(coef_H))
}




omega_t <- function(loadings, min_loading = .2) {
  
  loadings_group <- loadings[, -1, drop=FALSE]
  
  theta <- 1 - rowSums(loadings^2)  # model-implied unique variances
  
  # for the general factor
  outp <- sum(colSums(loadings)^2) / (sum(colSums(loadings)^2) + sum(theta))
  
  # for the group factors
  for (lupe in 1:ncol(loadings_group)) {
    
    group_TF <- abs(loadings_group[,lupe]) >= min_loading
    
    dum <- cbind( loadings[group_TF,1], loadings_group[group_TF,lupe])
    
    # 2023 Dueber, Toland - A Bifactor Approach to Subscore Assessment  p 224 formula 4
    outp <- c(outp, 
              ( sum(dum[,1])^2 + sum(dum[,2])^2 ) / 
                ( sum(dum[,1])^2 + sum(dum[,2])^2 + sum((theta*group_TF)) ) )
  }
  return(invisible(outp))
}



omega_h <- function(loadings, min_loading = .2) {
  
  loadings_group <- loadings[, -1, drop=FALSE]
  
  theta <- 1 - rowSums(loadings^2)  # model-implied unique variances
  
  # for the general factor
  outp <- sum(loadings[,1])^2 / (sum(colSums(loadings)^2) + sum(theta))
  
  # for the group factors
  for (lupe in 1:ncol(loadings_group)) {
    
    group_TF <- abs(loadings_group[,lupe]) >= min_loading
    
    dum <- cbind( loadings[group_TF,1], loadings_group[group_TF,lupe])
    
    # 2023 Dueber, Toland - A Bifactor Approach to Subscore Assessment  p 224 formula 5
    outp <- c(outp, 
              ( sum(dum[,2])^2 / 
                  ( sum(dum[,1])^2 + sum(dum[,2])^2 + sum(theta*group_TF)) ) )
  }
  return(invisible(outp))
}





# AI Overview for  r code for Schmid-Leiman Iterative Bi-factor Difference- based Target Rotation Algorithm
# 
# The Schmid-Leiman Iterative Bi-factor Difference-based Target Rotation Algorithm (SLiD) 
# is a state-of-the-art exploratory factor analysis method used to establish an 
# empirical, factor-specific cutoff for defining target matrices based on loading 
# differences.Because this algorithm relies on evaluating the largest drop 
# (one-lagged difference) in sorted squared loadings, you can implement the 
# full SLiD procedure in R using the base framework, psych package, and GPArotation 
# package. Complete R Implementation of SLiDCopy and execute this complete, 
# self-contained R script. It defines the SLiD algorithm function and applies it to 
# a sample datase.

# Algorithm Breakdown
# 
# Initial Alignment: The function utilizes psych::schmid to generate the foundational 
# hierarchical structure from first-order oblique factors.
# 
# Empirical Cutoff Rule: Instead of enforcing arbitrary cross-loading thresholds 
# (e.g., 0.10 or 0.20), SLiD dynamically evaluates individual factor vectors. 
# It maps out the sharpest mathematical drop in sorted squared loadings to 
# differentiate true indicator items from background noise.
# 
# Partially Specified Target Rotation: Using the computed GPArotation::pstQ system, 
# non-salient paths are held to zero while dominant indicators and cross-loadings 
# are left unconstrained (NA), resulting in a cleaner bi-factor structure.
# 
# To customize the script for your workflow, let me know:
#   Are your factors expected to be orthogonal or oblique?
#   What is the sample size and number of items in your dataset?
#   Are you planning to compare this against standard Direct Schmid-Leiman (DSL) or CFA implementations?

# ?? Google AI got the code from 
# "2019 - Improving Bi-Factor Exploratory Modeling Empirical Target Rotation Based on Loading Differences" ??


# -------------------------------------------------------------------------
# R Implementation: Schmid-Leiman Iterative Difference-Based Target Rotation (SLiD)
# -------------------------------------------------------------------------

#' SLiD: Schmid-Leiman Iterative Bi-factor Difference-based Target Rotation
#' 
#' @param data A matrix or dataframe of observed variables, or a correlation matrix.
#' @param n_factors Number of group/specific factors to extract.
#' @param max_iter Maximum number of outer target matrix update iterations.
#' @param tol Convergence tolerance (change in target matrix elements).
#'
slid_rotation <- function(data, n_factors = 3, max_iter = 100, tol = 1e-5,
                          fm = 'minres', rotate='oblimin') {
  
  # Ensure we have a correlation matrix
  if (is.matrix(data) && all(diag(data) == 1) && isSymmetric(data)) {
    R <- data
  } else {
    R <- cor(data, use = "pairwise.complete.obs")
  }
  
  # Step 1: Compute initial baseline Schmid-Leiman solution
  # This serves as the starting point for group factors
  sl_init <- psych::schmid(R, nfactors = n_factors, 
                           fm = fm, rotate=rotate, plot=FALSE)
  
  
  loadings_for_BIF <- sl_init$oblique[,1:n_factors]
  
  
  # Extract initial group factor loadings (exclude general factor column)
  # psych::schmid returns general factor in column 1, group factors in subsequent columns
  group_loadings <- sl_init$sl[, 2:(n_factors + 1)]
  p <- nrow(group_loadings)
  
  # Initialize the Target matrix (NA means unrestricted/free parameter)
  target_matrix <- matrix(NA, nrow = p, ncol = n_factors)
  old_target <- target_matrix
  
  # Main SLiD Iterative Loop
  for (iter in 1:max_iter) {
    
    # Calculate empirical cutoffs for each factor based on one-lagged differences
    for (j in 1:n_factors) {
      loadings_j <- group_loadings[, j]
      sq_loadings <- loadings_j^2
      
      # Sort squared loadings in descending order
      sorted_idx <- order(sq_loadings, decreasing = TRUE)
      sorted_sq <- sq_loadings[sorted_idx]
      
      # Compute first prominent one-lagged difference
      differences <- -diff(sorted_sq)
      
      if (length(differences) > 0) {
        max_diff_idx <- which.max(differences)
        # Cutoff is the average threshold separating major and minor loads
        cutoff_val <- sqrt(sorted_sq[max_diff_idx + 1])
      } else {
        cutoff_val <- 0.10 # Fallback default
      }
      
      # Construct target entries for factor j
      # Values below the empirical difference cutoff are targeted to 0
      target_matrix[, j] <- ifelse(abs(loadings_j) <= cutoff_val, 0, NA)
    }
    
    # Handle absolute strictness check for convergence
    if (iter > 1 && all(is.na(target_matrix) == is.na(old_target)) && 
        all(target_matrix[!is.na(target_matrix)] == old_target[!is.na(old_target)])) {
      # cat("SLiD Converged successfully at iteration:", iter, "\n")
      break
    }
    
    old_target <- target_matrix
    
    # Apply Partially Specified Procrustes/Target Rotation (pstQ or pstT) via GPArotation
    # We rotate the group factor space conditionally based on our target matrix
    tryCatch({
      rotated_solution <- GPArotation::pstQ(group_loadings, target = target_matrix)
      group_loadings <- rotated_solution$loadings
    }, error = function(e) {
      # warning("Target rotation matrix issue encountered during iteration ", iter, ". Using last stable state.")
    })
  }
  
  # Construct final Bifactor Matrix output layout
  # Re-attach general factor from baseline for complete structural representation
  final_bifactor_matrix <- cbind(sl_init$sl[, 1, drop = FALSE], group_loadings)
  colnames(final_bifactor_matrix) <- c("General", paste0("Group_", 1:n_factors))
  
  return(list(
    bifactor_loadings = final_bifactor_matrix,
    final_target = target_matrix,
    iterations = iter,
    loadings_for_BIF = loadings_for_BIF
  ))
}





