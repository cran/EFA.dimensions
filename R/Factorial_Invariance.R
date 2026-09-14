


Factorial_Invariance <- function(data, group, 
                                 LV_model=NULL, 
                                 LV_keys=NULL, 
                                 LV_names=NULL, 
                                 LV_resid_correls=NULL, 
                                 estimator = 'ML', 
                                 ordered = FALSE, 
                                 disp_options = list(models = c('Configural', 'Scalar'),
                                                     coefs = c('std', 'LV_means')),  
                                 verbose=TRUE) {
  
  ###############################  data  setup   ####################################
  
  data <- MISSING_DROP(data)
  
  varnames <- colnames(data)[colnames(data) != group]
  

  #############################  argument checks   ##################################

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
  
  if (!(estimator %in% c('ML', 'GLS', 'WLS', 'DWLS', 'ULS', 'DLS', 'PML', 'MLM', 
                         'MLMVS', 'MLMV', 'MLF', 'MLR', 'WLSM', 'WLSMVS', 'WLSMV', 
                         'ULSM', 'ULSMVS', 'ULSMV')))
    stop('The provided estimator, ', estimator, ', is not one of the possible estimators.')
  

  ##################### check model / create lavaan model syntax   ##################
  
   lav_model <- lavaan_model(varnames = varnames, model=LV_model, keys=LV_keys, 
                            LV_names=LV_names, resid_correls=LV_resid_correls)
  
  
  ###################################  lavaan   #####################################
  
  fit_coef_prefs <- c('rmsea','srmr','cfi','tli','aic','bic')
  
  # npar fmin chisq df pvalue baseline.chisq baseline.df baseline.pvalue 
  # cfi tli nnfi rfi nfi pnfi ifi rni logl unrestricted.logl aic bic ntotal 
  # bic2 rmsea rmsea.ci.lower rmsea.ci.upper rmsea.ci.level rmsea.pvalue 
  # rmsea.close.h0 rmsea.notclose.pvalue rmsea.notclose.h0 rmr rmr_nomean
  # srmr srmr_bentler srmr_bentler_nomean crmr crmr_nomean 
  # srmr_mplus srmr_mplus_nomean cn_05 cn_01 gfi agfi pgfi mfi ecvi
  
  # get fit coefs for each group separately
  group_values <- unique(data[,group])
  fits_by_group <- c()
  for (lupe in 1:length(group_values)) {
    
    mod_all <- cfa(model = lav_model, 
                   data = subset(data, eval(parse(text=group)) == group_values[lupe]),
                   estimator = estimator, ordered = ordered)
    
    fits_by_group <- rbind(fits_by_group, fitMeasures(mod_all)[fit_coef_prefs])
  }
  rownames(fits_by_group) <- group_values
  
  
  mod_Configural <- cfa(model = lav_model, data=data, group=group, 
                        estimator = estimator, ordered = ordered)
  
  mod_Metric <- cfa(model = lav_model, data=data, group=group,
                    group.equal = c('loadings'), 
                    estimator = estimator, ordered = ordered)
  
  mod_Scalar <- cfa(model = lav_model, data=data, group=group,
                    group.equal = c('loadings','intercepts'), 
                    estimator = estimator, ordered = ordered)
  
  mod_Strict <- cfa(model = lav_model, data=data, group=group,
                    group.equal = c('loadings','intercepts', 'residuals'),
                    estimator = estimator, ordered = ordered)
  
  mod_LV_vars <- cfa(model = lav_model, data=data, group=group,
                     group.equal = c('loadings','intercepts', 'residuals', 
                                     'lv.variances'),
                     estimator = estimator, ordered = ordered)
  
  LV_noms <- lavNames(mod_Configural, type = "lv")
  
  if (length(LV_noms) > 1) {
    mod_LV_covars <- cfa(model = lav_model, data=data, group=group,
                         group.equal = c('loadings','intercepts', 'residuals', 
                                         'lv.variances', 'lv.covariances'),
                         estimator = estimator, ordered = ordered)
    
    mod_LV_means <- cfa(model = lav_model, data=data, group=group,
                        group.equal = c('loadings','intercepts', 'residuals', 
                                        'lv.variances', 'lv.covariances', 'means'),
                        estimator = estimator, ordered = ordered)
  } else {
    mod_LV_covars <- NULL 
    
    mod_LV_means <- cfa(model = lav_model, data=data, group=group,
                        group.equal = c('loadings','intercepts', 'residuals', 
                                        'lv.variances', 'means'),
                        estimator = estimator, ordered = ordered)
  } 
  
  # # measurementInvariance doesn't do equal factor variance. But, this can be accomplished as follows
  # fit.ef <- cfa(model.mdd, data = Data, group = "sex", meanstructure = TRUE) # equal form
  # fit.efl <- update(fit.ef, group.equal = c("loadings")) # equal factor laodings
  # fit.eii <- update(fit.efl, group.equal = c("loadings", "intercepts")) # equal indicator intercepts
  # fit.eir <- update(fit.eii, group.equal = c("loadings", "intercepts", "residuals")) # equal indicator error variances
  # fit.fv <- update(fit.eir, group.equal = c("loadings", "intercepts", "residuals", "lv.variances")) # equal factor variances
  # fit.fm <- update(fit.fv, group.equal = c("loadings", "intercepts", "residuals", "lv.variances", "means")) # equal latent means
  
  # std_sol <- standardizedSolution(mod_LV_means)
  # 
  # # Filter for factor loadings (op == "=~")
  # std_loadings <- std_sol[std_sol$op == "=~", ]
  # print(std_loadings)
  
  
  # loadings_std <- lavInspect(mod_Strict, what = "std")$lambda
  # cat('\n\nStandardized loadings\n\n'); print(round(loadings_std,2), print.gap=4)
  # 
  # loadings_std <- lavInspect(mod_LV_means, what = "std")$lambda
  # cat('\n\nStandardized loadings\n\n'); print(round(loadings_std,2), print.gap=4)
  # 
  
  
  # there is/was a problem running summary on an semTools::compareFit object, as described
  
  # https://groups.google.com/g/lavaan/c/RqKf4_wQ1e4
  # 
  # https://github.com/simsem/semTools/issues/116
  # 
  # https://stackoverflow.com/questions/63065696/error-in-getmethodsummary-signature-fitdiff/63600180#63600180
  
  # it works ok outside of a function, but not within my function
  # the problem is with the object type, S4 vs S3
  
  # this code extracts the semTools::compareFit object into a useable list
  # it should appear before the semTools::compareFit command
  # NOTE: it does not work when just running the commands by themselves, outside of 
  #       package building
  summary.FitDiff <- function(object){
    output <- list(fit = object@fit, fit.diff = object@fit.diff,
                   chisq_diffs = object@nested)
    return(invisible(output))
  }
  
  
  if (length(lavNames(mod_Configural, type = "lv")) > 1) {
    CompareAll <- semTools::compareFit(mod_Configural, 
                                       mod_Metric,
                                       mod_Scalar,
                                       mod_Strict,
                                       mod_LV_vars,
                                       mod_LV_covars,
                                       mod_LV_means)
  } else {
    CompareAll <- semTools::compareFit(mod_Configural,  
                                       mod_Metric,
                                       mod_Scalar,
                                       mod_Strict,
                                       mod_LV_vars,
                                       mod_LV_means)
  }
  
  # # to prevent display of output (re: from summary(CompareAll ...)
  ## not needed when using the above summary.FitDiff function
  # quiet_sum <- function(x) {sink(tempfile()); on.exit(sink()); invisible(summary(x))}
  # 
  # mod_comps <- quiet_sum(summary(CompareAll))
  
  mod_comps <- summary(CompareAll)
  
  # model_fits <- mod_comps@fit
  model_fits <- mod_comps$fit
  model_fits <- model_fits[,fit_coef_prefs]
  
  # model_fit_diffs <- mod_comps@fit.diff
  model_fit_diffs <- mod_comps$fit.diff
  model_fit_diffs <- model_fit_diffs[,fit_coef_prefs] 
  
  # chisq_diffs <- mod_comps@nested
  chisq_diffs <- mod_comps$chisq_diffs
  chisq_diffs <- chisq_diffs[,c('Df','AIC','BIC','Chisq','Chisq diff','RMSEA',
                                'Df diff','Pr(>Chisq)')]
  
  
  ###################################  LV scores  ###################################
  
  LV_scores_outp <- LV_grp_means_outp <- LV_grp_sds_outp <- LV_grp_Ns_outp <- NULL
  
  if (!'none' %in% disp_options$models & 'LV_means' %in% disp_options$coefs) {
    
    LV_noms <- lavNames(mod_all, type = "lv")
    
    N_LVs <- length(LV_noms)
    
    group_noms <- group_values  # lavInspect(mod_Scalar, "group.label")
    
    N_groups <- length(group_noms)
    
    for (lupe in disp_options$models) {
      
      obj <- eval(parse(text = paste('mod_', lupe, sep='')))
      
      # extract latent variable scores
      LV_scores_list <- lavPredict(obj)
      
      LV_grp_means <- LV_grp_sds <- LV_grp_Ns <- matrix(NA, N_LVs, N_groups)
      
      LV_scores <- data.frame(matrix(ncol = (N_LVs + 1), nrow = 0))
      colnames(LV_scores) <- c(LV_noms, group)
      
      for (lupe_grps in 1:N_groups) {
        
        don <- data.frame(LV_scores_list[[lupe_grps]])
        
        don[group] <- group_noms[lupe_grps]
        
        LV_scores <- rbind(LV_scores, don)
        
        LV_grp_means[,lupe_grps] <- colMeans(LV_scores_list[[lupe_grps]])
        
        LV_grp_sds[,lupe_grps] <- sapply(data.frame(LV_scores_list[[lupe_grps]]), sd, na.rm = TRUE)
        
        LV_grp_Ns[,lupe_grps] <- nrow(don)
      }
      
      rownames(LV_grp_means) <- rownames(LV_grp_sds) <- rownames(LV_grp_Ns) <- LV_noms
      colnames(LV_grp_means) <- colnames(LV_grp_sds) <- colnames(LV_grp_Ns) <- group_noms
      
      LV_scores_outp[[lupe]] <- LV_scores
      LV_grp_means_outp[[lupe]] <- LV_grp_means
      LV_grp_sds_outp[[lupe]] <- LV_grp_sds
      LV_grp_Ns_outp[[lupe]] <- LV_grp_Ns
    }
  }
  
  
  ###################################  output  ######################################
  
  if (verbose) {
    
    cat('\n\nFactorial_Invariance\n')
    
    cat('\nThe model:\n')
    cat(lav_model, '\n')
    
    cat('\n\nThe levels of "group" and their Ns:')
    print(table(data[,group]))
    
    cat('\n\nFit coefficients for each group separately:\n')
    print(round(fits_by_group, 3), print.gap=4)
    
    cat('\n\nInvariance Model Fit Coefficients:\n')
    print(round(model_fits, 3), print.gap=4)
    
    cat('\n\nDifferences in Invariance Model Fit Coefficients:\n')
    print(round(model_fit_diffs, 3), print.gap=4)
    
    cat('\n\nChi-Squared Difference Tests:\n')
    print(round(chisq_diffs, 3), print.gap=4)
    
    if (!'none' %in% disp_options$models) {
      
      for (lupe in disp_options$models) {
        
        cat('\n\nCoefficients for', lupe)
        
        obj <- eval(parse(text = paste('mod_', lupe, sep='')))
        
        if ('raw' %in% disp_options$coefs) {
          
          cat('\n\nUntandardized loadings for each group (', lupe, '):\n\n', sep='')
          
          print(lapply(lapply(lavInspect(obj, what = "est"), function(x) x$lambda), round, digits = 2), print.gap=4)
        }
        
        if ('std' %in% disp_options$coefs) {
          
          cat('Standardized loadings for each group (', lupe, '):\n\n', sep='')
          
          print(lapply(lapply(lavInspect(obj, what = "std"), function(x) x$lambda), round, digits = 2), print.gap=4)
        }
        
        if ('phi' %in% disp_options$coefs) {  # psi = latent variable covariance matrix.
          
          cat('Latent variable intercorrelations for each group (', lupe, '):\n\n', sep='')
          
          print(lapply(lavInspect(obj, what = "cor.lv"), round, digits = 2), print.gap=4)
        }
        
        if ('LV_means' %in% disp_options$coefs) {
          
          # warning if the latent means are the same for the groups 
          # (e.g., when intercepts & loadings were not constricted to be equal)
          # need at least Scalar, & not Metric or Configural
          if (all(round(LV_grp_means_outp[[lupe]],4) == 0)) {
            cat('\n\nWARNING: The latent trait means for the groups are apparently the same.')
            cat('\nThis occurs in a model where the intercepts & loadings were not constrained to be equal.')
            cat('\nConsider using a "Scalar" model instead.\n')
          }
          
          cat('Latent variable group means (', lupe, '):\n\n', sep='')
          print(round(LV_grp_means_outp[[lupe]],2), print.gap=4)
          
          cat('\nLatent variable group standard deviations (', lupe, '):\n\n', sep='')
          print(round(LV_grp_sds_outp[[lupe]],2), print.gap=4)
          
          cat('\nLatent variable group Ns (', lupe, '):\n\n', sep='')
          print(LV_grp_Ns_outp[[lupe]], print.gap=4)
          
          # library(DFA.CANCOR)
          # ?GROUP.PROFILES
          # 
          # GROUP.DIFFS(LV_scores, GROUPS = 'group', DV = LV_noms[1], var.equal=FALSE, p.adjust.method="fdr")
          # 
          # GROUP.DIFFS(LV_scores, GROUPS = 'group', DV = 'pos', var.equal=FALSE, p.adjust.method="fdr")
        }
      }
    }
  }
  
  output <- list(lavaan_model_syntax = lav_model,
                 fits_by_group=fits_by_group, 
                 mod_Configural = mod_Configural,
                 mod_Metric     = mod_Metric,
                 mod_Scalar     = mod_Scalar,
                 mod_Strict     = mod_Strict,
                 mod_LV_vars    = mod_LV_vars,
                 mod_LV_covars  = mod_LV_covars,
                 mod_LV_means   = mod_LV_means,
                 model_fits=model_fits, 
                 model_fit_diffs=model_fit_diffs, 
                 chisq_diffs=chisq_diffs,
                 LV_scores = LV_scores)
  
  class(output) <- "Factorial_Invariance"
  
  return(invisible(output))  
}

