\name{PLOT_Invariance}
\alias{PLOT_Invariance}
\title{Plots of multiple groups factor analysis invariance statistics}
\description{Plots factor loadings, slopes and intercepts, and/or latent variable
distributions from multiple groups confirmatory factor analysis (CFA) models
and exploratory structural equation models (ESEM)}
\usage{PLOT_Invariance(model_object,
                invar_model = 'Scalar',
                plot_types = c('loadings', 'ints_slopes', 'LV_distribs'),
                loadings_std = TRUE, loadings_by_factor = FALSE)}
\arguments{
  \item{model_object}{
  \code{}A "Factorial_Invariance" function output object, or a lavaan 
  multiple groups cfa output object.}

  \item{invar_model}{
  \code{} The kind of invariance model for the plots when model_object is from 
  the Factorial_Invariance function. The options are: 'Configural', 
  'Metric', 'Scalar', 'Strict', 'LV_vars', 'LV_covars', and 'LV_means'.
  (This argument can be ignored when model_object is not from the
  Factorial_Invariance function.)}

  \item{plot_types}{
  \code{}The desired kind(s) of plots. The possibilities are
  \cr (1) 'loadings', for factor loadings, 
  \cr (2) 'ints_slopes', for intercepts and slopes, and/or 
  \cr (3) 'LV_distribs', for latent variable distributions.}
  	    
  \item{loadings_std}{
  \code{}Should standardized loadings be used for 'loadings' plots? TRUE or FALSE}
  
  \item{loadings_by_factor}{
  \code{}For 'loadings' plots, should a separate plot be produced for each 
  factor (TRUE)? If FALSE, then the loadings for different factors will 
  appear in the same plots.}
}
\value{A list containing the following components (if plotted using this function):
loadings, estimates, & LV_scores.
}
\author{Brian P. O'Connor }
\examples{
RSE_output <-
Factorial_Invariance(data = data_RSE_sex, group = 'gender', 
                     LV_keys = c(1, 1, 2, 1, 2, 1, 1, 2, 2, 2), 
                     LV_names = c('pos', 'neg'))

PLOT_Invariance(model_object = RSE_output, 
                invar_model = 'Configural',
                plot_types = c('loadings', 'ints_slopes'),
                loadings_std = TRUE, loadings_by_factor = TRUE) 
\donttest{
HS_1939_output <-
  Factorial_Invariance(data = data_HS_1939, group = 'school',
                       LV_keys = c(1,1,1, 2,2,2, 3,3,3),
                       LV_names = c('visual', 'textual', 'speed'))

PLOT_Invariance(model_object = HS_1939_output, 
                invar_model = 'Scalar',
                plot_types = c('LV_distribs')) 
}
}
