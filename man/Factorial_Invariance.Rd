\name{Factorial_Invariance}
\alias{Factorial_Invariance}
\title{Factor Model Invariance}
\description{Assesses the equivalence of measurement and structural parameters 
of a CFA or ESEM factor model across multiple groups}
\usage{
Factorial_Invariance(data, 
                     group,
                     LV_model=NULL,
                     LV_keys=NULL,
                     LV_names=NULL, 
                     LV_resid_correls=NULL,
                     estimator = 'ML',
                     ordered = FALSE,
                     disp_options = list(models = c('Configural', 'Scalar'),
                                         coefs = c('std', 'LV_means')),
                     verbose=TRUE)
}
\arguments{
  \item{data}{
  \code{}A dataframe with the model variables and the group factor.
  \cr Example: data = data_RSE_sex}

  \item{group}{
  \code{}The name of the group variable in data.
  \cr Example: group = 'gender'}

  \item{LV_model}{
  \code{}(optional) lavaan package CFA model syntax. Example:
  \cr model = '
  \cr pos_items =~ Q1 + Q2 + Q4 + Q6 + Q7
  \cr neg_items =~ Q3_R + Q5_R + Q8_R + Q9_R + Q10_R'}

  \item{LV_keys}{
  \code{}(optional) A vector or dataframe indicating which item (i.e., which variable
  in data) goes with which latent variable (which can be named or numbered).
  LV_keys is required when LV_model is not provided. If LV_keys is a dataframe,
  then the variable names must be in the first column and the factor memberships
  must be in the second column. Example of a vector entry:
  \cr LV_keys = c(1, 1, 2, 1, 2, 1, 1, 2, 2, 2)}

  \item{LV_names}{
  \code{}(optional) A vector with the preferred names for the latent variables. 
  \cr Example: LV_names = c('pos', 'neg')}

  \item{LV_resid_correls}{
  \code{}(optional) A vector with the pairs of correlated error terms (if any). 
  \cr Example: LV_resid_correls = c('Q1 ~~ Q2', 'Q9 ~~ Q10')}

 \item{estimator}{
  \code{}The name of the lavaan estimator to be used in the analyses.
  \cr\cr The options for basic estimators (continuous data) are 
  ML, GLS, WLS, DWLS, ULS, DLS, and PML.
  \cr\cr The options for robust estimators are MLM, MLMVS, MLMV, MLF, MLR, 
  WLSM, WLSMVS, WLSMV, ULSM, ULSMVS, and ULSMV.
  \cr\cr Example: estimator = 'ML'}

  \item{ordered}{
  \code{}(optional) are the variables in data (if provided) ordered? 
    For example, Likert scale responses are usually ordered. 
    The possible values for ordered are TRUE, FALSE, or a vector with 
    the names of the ordered variables.
}

  \item{disp_options}{
  \code{}(optional) A list with display (print) options for model coefficients. 
  The list elements are 'models' and 'coefs'. 
  \cr\cr The possible entries for 'models' are: 'mod_Configural', 'mod_Metric', 
  'mod_Scalar', 'mod_Strict', 'mod_LV_vars', 'mod_LV_covars', & 'mod_LV_means'.
  \cr\cr The possible entries for 'coefs' are: 'raw', 'std', 'phi', & 'LV_means'
  \cr\cr See the Details section below for information on each element.}
  
  \item{verbose}{
  \code{}(optional) Should detailed results be displayed in console? 
  \cr TRUE (default) or FALSE}
}
\details{
  This function assesses the equivalence of the measurement and structural 
  parameters of a confirmatory factor analysis (CFA) model, or of an exploratory
  structural equation model (ESEM), across multiple groups.
  It is a wrapper function that uses functions from the lavaan (Rosseel, 2012) 
  and semTools (Jorgensen et al., 2026) packages.
  
  The measurement model pertains to the measurement characteristics of 
  the indicators (observed measures) and thus consists of the factor 
  loadings, intercepts, and residual variances. The evaluations of the 
  equivalence of these parameters across groups are tests of measurement 
  invariance (Brown, 2015).

  Tests of the structural parameters involve evaluations of the latent 
  variables themselves and focus on the factor variances, covariances, 
  and latent means. These parameters are (estimated) features of the 
  population from which a sample was drawn. The evaluations of the 
  equivalence of these parameters across groups are tests of population 
  heterogeneity: Does the dispersion, the interrelationships, and the 
  levels of the factors vary across groups?

  Brown (2015) stated that it is more prudent for model evaluation to 
  work upward from the least restricted solution (equal form) to 
  determine if further tests of measurement invariance and population 
  heterogeneity are warranted. Here is the recommended sequence of 
  multiple groups invariance evaluations that is used in the present 
  function, with the function model names in bold:

	\itemize{
		
    \item\strong{Configural} invariance, or "equal form" invariance, 
    tests whether the number of factors and the pattern of 
    indicator factor loadings are identical across groups.

    \item\strong{Metric} invariance, or "equal factor loadings" invariance, 
    tests the equality of the factor loadings across groups.

    \item\strong{Scalar} invariance, "equal intercepts" invariance, or 
    "strong factorial" invariance, tests whether the indicator intercepts 
    are equal across groups.

    \item\strong{Strict} invariance, or "equal residual" invariance, 
    tests whether the indicator residual variances are equal across groups.

    \item\strong{LV_vars} tests whether the factor variances are equal across groups.

    \item\strong{LV_covars} tests whether the factor covariances are equal across groups.

    \item\strong{LV_means} tests whether the latent variable means are equal across groups.
    }
}
\value{A list containing the following components:
\item{lavaan_model_syntax}{The lavaan CFA model syntax}
\item{fits_by_group}{The model fit coefficients for each group separately}
\item{mod_Configural}{The lavaan cfa output object for the configural invariance model}
\item{mod_Metric}{The lavaan cfa output object for the metric invariance model}
\item{mod_Scalar}{The lavaan cfa output object for the scalar invariance model}
\item{mod_Strict}{The lavaan cfa output object for the strict invariance model}
\item{mod_LV_vars}{The lavaan cfa output object for the model that restricts the 
latent variable variances to be equal}
\item{mod_LV_covars}{The lavaan cfa output object for the model that restricts 
the latent variable covariances to be equal}
\item{mod_LV_means}{The lavaan cfa output object for the model that restricts the 
latent variable means to be equal}
\item{model_fits}{The invariance model fit coefficients}
\item{model_fit_diffs}{Differences in the invariance model fit coefficients}
\item{chisq_diffs}{Chi-Squared model differences tests}
\item{LV_scores}{A list with scores on the latent variables for each of the models in disp_options}
}
\references{       
       {Brown, T. A. (2015). \emph{Confirmatory Factor Analysis for Applied 
       Research (2nd ed.)}. New York: The Guilford Press.}
\cr\cr {Jorgensen, T. D., Pornprasertmanit, S., Schoemann, A. M., & 
       Rosseel, Y. (2026). semTools: Useful tools for structural equation 
       modeling. R package version 0.5-8. 
       Retrieved from https://CRAN.R-project.org/package=semTools}
\cr\cr {Milfont, T. L., & Fischer, R. (2015). Testing measurement invariance across 
        groups: Applications in cross-cultural research. 
        \emph{International Journal of Psychological Research, 3,} 111130.}
\cr\cr {Rosseel, Y. (2012). lavaan: An R package for structural equation modeling. 
        \emph{Journal of Statistical Software, 48,} 136.}
\cr\cr {Stark, S., Chernyshenko, O. S., & Drasgow, F. (2006). Detecting differential 
        item functioning with confirmatory factor analysis and item response theory: 
        Toward a unified strategy. \emph{Journal of Applied Psychology, 91(6),} 12921306.}
\cr\cr {Stein, J. A., Lee, J. W., & Jones, P. S. (2006). Assessing cross-cultural 
        differences through use of multiple-group invariance analyses. 
        \emph{Journal of Personality Assessment, 87(3),} 249258.}
\cr\cr {Tan, T. (2024) Frequentist and Bayesian factorial invariance using R. 
        \emph{Practical Assessment, Research, and Evaluation. 29(8),} 1-34.}
}
\author{Brian P. O'Connor}
\examples{
model_RSE <- '
   pos_items =~ Q1 + Q2 + Q4 + Q6 + Q7
   neg_items =~ Q3_R + Q5_R + Q8_R + Q9_R + Q10_R '

# using LV_model
Factorial_Invariance(data = data_RSE_sex, group = 'gender', LV_model = model_RSE)
\donttest{
# using LV_keys
Factorial_Invariance(data = data_RSE_sex, group = 'gender', 
                     LV_keys = c(1, 1, 2, 1, 2, 1, 1, 2, 2, 2), 
                     LV_names = c('pos', 'neg'))

# another way of using LV_keys
Factorial_Invariance(data = data_RSE_sex, group = 'gender', 
                     LV_keys = c( Q1 = 1, Q2 = 1, Q3_R = 2, Q4 = 1,   Q5_R = 2, 
                                  Q6 = 1, Q7 = 1, Q8_R = 2, Q9_R = 2, Q10_R = 2) )

# for an ESEM model
# first, use the ESEM function to obtain the esem_model_syntax (without the group variable, gender)
esem_output <- ESEM(data = subset(data_RSE_sex, select = -c(gender)), Nfactors = 2)

Factorial_Invariance(data = data_RSE_sex, group = 'gender',
                     LV_model = esem_output$esem_model_syntax)
        
         
# See this web page for additional Examples:

# https://oconnor-psych.ok.ubc.ca/vignettes/EXAMPLES_vignettes.html
}
}
