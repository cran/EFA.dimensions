\name{OMEGA}
\alias{OMEGA}
\title{Omega internal consistency reliability coefficients}
\description{Total and hierarchical omega internal consistency reliability 
coefficients computed using multiple possible methods}
\usage{OMEGA(data, corkind = 'pearson',  
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
      min_loading = .2, display = 2)} 
\arguments{
  \item{data}{
  \code{}An all-numeric dataframe where the rows are cases & the columns are 
  the variables, or a correlation matrix with ones on the diagonal. The function 
  internally determines whether the data are a correlation matrix.}

  \item{corkind}{
  \code{}The kind of correlation matrix to be used if data is not a correlation matrix. 
  The options are 'pearson', 'kendall', 'spearman', 'gamma', and 'polychoric'.}

  \item{Nfactors}{
  \code{}(optional) The number of factors to extract, including the general factor.
  \cr Nfactors is = 1 (for the general factor) + the number of group factors.}

  \item{bifactor_kind}{
  \code{}The bifactor method(s) to be used. The options 
  are 'SL', 'SLiD', 'DSL', 'bifactorQ', 'bifactorT', 
  'bigeominQ', 'bigeominT', 'CFA', and 'ESEM'. Multiple methods can be specified.}

  \item{EFA_options}{
  \code{}(optional) A list with EFA options when rawdata is provided. The list elements
  include values for 'extraction' and 'rotation'.
  \cr\cr The possibilities for extraction are: 'minres'(the default), 'alpha', 'fullinfo', 
  'gls', 'image', 'ml', 'ols', 'paf', 'uls', and 'wls'.
  \cr\cr The possibilities for rotation are:  
  'varimax', 'bentlerT', 'entropy', 'equamax', 'geominT', 'quartimax', 
  'promax', 'bentlerQ', 'geominQ', 'oblimin' (the default, 'oblimax', 'quartimin', 
  'simplimax', and 'none'.
  }
 
 \item{LV_options}{
  \code{}(optional) A list with cfa options when bifactor_kind = 'CFA'. 
  The possible list elements are 'group_keys',
  'estimator', 'rotation', 'resid_correls', 'LV_names', and 'ordered'. 
  See the Details section below for information on each element.}
  
  \item{schmid_options}{
  \code{}(optional) A list with schmid function (from the psych package) options 
  when bifactor_kind is one of 'SL', 'SLiD', or 'DSL'. The list elements
  include values for 'extraction' and 'rotation'.
  \cr\cr The possibilities for extraction are: 'minres'(the default), 'paf', 'pc', and 'ml'.
  \cr\cr The possibilities for rotation are: 'oblimin' (the default), 'simplimax', 
  'Promax', 'promax' and 'none'.}

  \item{GPA_options}{ 
  \code{}(optional) A list with options for when a GPArotation is used.
  See the Details section below.}
  
 \item{min_loading}{ 
  \code{} The minimum value of a group factor loading for an item to be
  considered to have a non-negligible contribution to a group factor.
  min_loading only plays a role in the computations for PUC and for
  some group factor statistics.}
                  
  \item{display}{
  \code{}The results to be displayed in the console: 
  \cr 0 = nothing
  \cr 1 = only the omega and scale coefficients
  \cr 2 (default) = detailed output, including the bifactor model loadings and statistics}
}
\details{

For the \strong{bifactor_kind} argument:

	\itemize{
		
    \item\strong{bifactorT} is an orthogonal bifactor rotation 
    designed for situations where a single, overarching global dimension is expected
    alongside separate sub-domains, and when all factors should be 
    uncorrelated.
        
    \item\strong{bifactorQ} is an oblique bifactor rotation designed for when 
    a strong, overarching global dimension is expected alongside separate 
    sub-domains that are allowed to correlate with each other.

    \item\strong{bigeominT} is an orthogonal bifactor rotation 
    designed for situations where a single, overarching global dimension is expected
    alongside separate sub-domains, and when all factors should be 
    uncorrelated.
        
    \item\strong{bigeominQ} is an oblique bifactor rotation designed for when 
    a strong, overarching global dimension is expected alongside separate 
    sub-domains that are allowed to correlate with each other.

    \item\strong{SL} is an oblique bifactor rotation designed for 
    where there is a broad, overarching factor 
    alongside sub-domains that are allowed to overlap. 
    
    \item\strong{SLiD} an orthogonal bifactor rotation method designed for 
    exploratory bifactor analysis that tries to ensure that 
    once the general factor's variance is pulled out, each item loads 
    onto only one specific group factor. 
        
    \item\strong{DSL} entropy is an orthogonal factor rotation 
    that pushes factor loadings to be either 
    strongly dominant (close to 1.0) or cleanly absent (close to 0.0), 
    reducing the overall informational "noise" of the matrix.
    
    \item\strong{CFA} confirmatory bifactor analyses are conducted
    using the cfa function from the lavaan package. All items are
    assumed to load on the general factor. The items that load on
    the group factors are indicated using "group_keys" for the
    LV_options argument. 
    
    \item\strong{ESEM} rxploratory structural equation modeling bifactor 
    analyses conducted
    using the cfa function from the lavaan package. All items are
    assumed to load on the general factor. The items that load on
    the group factors can be indicated using "group_keys" for the
    LV_options argument. \cr
  }


For the \strong{LV_options} argument:

  When bifactor_kind = 'CFA' or 'ESEM', then LV_options is a list with information for 
  these analyses. The possible list elements are 'group_keys', 
  'estimator', 'rotation', 'resid_correls', 'LV_names', and 'ordered'.

	\itemize{
		
    \item\strong{group_keys} is a vector indicating which item (i.e., which variable
    in data) goes with which group factor (which can be named or numbered). \cr
    Example: group_keys = c(1, 1, 2, 1, 2, 1, 1, 2, 2, 2)
  
    \item\strong{estimator} is the name of the lavaan estimator to be used in the analyses.
    The options for basic estimators (continuous data) are 
    ML, GLS, WLS, DWLS, ULS, DLS, and PML.
    The options for robust estimators are MLM, MLMVS, MLMV, MLF, MLR, 
    WLSM, WLSMVS, WLSMV, ULSM, ULSMVS, and ULSMV. \cr
    Example: estimator = 'ML'

    \item\strong{rotation} is the kind of lavaan bifactor rotation. The two
    options are 'bigeomin' or 'biquartimin'.

    \item\strong{resid_correls} is a vector with the pairs of correlated error terms (if any). \cr 
    Example: resid_correls = c('Q1 ~~ Q2', 'Q9 ~~ Q10')
    
    \item\strong{LV_names} is a vector of names for the latent variables
     
    \item\strong{ordered} are the variables in rawdata (if provided) ordered? 
    For example, Likert scale responses are usually ordered. The function will
    produce notices when the data are not consistent with the "ordered" status.
  }

  Examples of LV_options:
  
    \preformatted{
  LV_options = list(group_keys = c(1, 1, 2, 1, 2, 1, 1, 2, 2, 2),
                     resid_correls = NULL, estimator = 'WLSMVS')
  }
    \preformatted{
  LV_options = list(group_keys = c(Q1=1, Q2=1, Q3_R=2, Q4=1,   Q5_R=2, 
                                    Q6=1, Q7=1, Q8_R=2, Q9_R=2, Q10_R=2)),
                     resid_correls = c('Q1 ~~ Q2', 'Q9 ~~ Q10'), 
                     estimator = 'WLSMVS')
}


  For the \strong{GPA_options} argument:

  A list with options for when a GPArotation is used. The list 
  elements are 'delta', 'epsilon', 'normalize', 'maxit', and 'randomStarts'.
  
  \itemize{
		
    \item\strong{delta} is a small positive number that is added to the squared  
    factor loadings in the (geomin) optimization denominator to prevent division by 
    zero and control the sharpness of the penalty near zero. 
    (Note: In MPLus, delta is called "epsilon".)
  
    \item\strong{epsilon} is the numeric rotational convergence tolerance threshold. 
    It controls how precise or strict the stopping rule is for the rotation algorithm.
    The gradient projection algorithm stops optimizing and assumes the 
    rotation has converged when the norm of the gradient falls below this value.

    \item\strong{normalize} specifies whether and how the loadings
    matrix is normalized prior to rotation and denormalized after rotation.
    The possibilities are FALSE (the default) for no normalization, or
    TRUE, for the Kaiser/Horst normalization.

    \item\strong{maxit} is the maximum number of iterations.

    \item\strong{randomStarts} is an integer that specifies the number of 
    random starting configurations used to optimize the factor rotation.
    It helps avoid getting trapped in local minima during factor rotation 
    optimization. It runs the gradient projection algorithm multiple timeseach 
    using a different random orthogonal or oblique starting rotation 
    matrix (Tmat) generated via QR decomposition.
    It returns the rotated factor loadings matrix and results corresponding 
    to the run with the lowest objective criterion value among all attempts.
  }


Run the following command for detailed descriptions of 
the omega-total and omega-hierarchical internal consistency reliability 
coefficients and other statistics produced by this function: \cr
\itemize{

  \item{RShowDoc("Coefficient_descriptions_vignettes", package = "EFA.dimensions")} \cr

  \item{vignette("Coefficient_descriptions_vignettes")} \cr
}
}
\value{A list with the omega coefficients, the factor loadings, and model fit statistics.
}
\references{ 
      {Bernaards, C. A., & Jennrich, R. I. (2005). Gradient Projection Algorithms 
       and Software for Arbitrary Rotation Criteria in Factor Analysis. 
       \emph{Educational and Psychological Measurement, 65(5),} 676-696.}
\cr\cr {Bernaards, C. A., & Jennrich, R. I. (2026). GPArotation: Gradient Projection 
       Factor Rotation. R package version 2026.4-1, 
       https://CRAN.R-project.org/package=GPArotation}
\cr\cr {Educational Content Team. (2026). McDonald's Omega. Cogn-IQ Encyclopedia. 
        https://pubscience.org/cqep.2025.0048}
\cr\cr {Flora, D. B. (2020). Your coefficient alpha is probably wrong, but which 
       	coefficient omega is right? A tutorial on using R to obtain better reliability estimates. 
       	\emph{Advances in Methods and Practices in Psychological Science, 3(4),} 484501.}
\cr\cr {Garcia-Garzon, E., Abad, F. J., & Garrido, L. E. (2021). On omega  
        hierarchical estimation: A Comparison of Exploratory Bi-Factor Analysis
        Algorithms. \emph{Multivariate Behavioral Research, 56(1),} 101-119.}
\cr\cr {Jennrich, R. I. (2018). Rotation. In P. Irwing, T. Booth, & D. J. 
        Hughes (Eds.), \emph{The Wiley handbook of psychometric testing: A 
        multidisciplinary reference on survey, scale and test development} (pp. 279304). 
 	      Wiley Blackwell. https://doi.org/10.1002/9781118489772.ch10}
\cr\cr {Kalkbrenner, M. T. (2024). Choosing Between Cronbachs Coefficient Alpha, 
        McDonalds Coefficient Omega, and Coefficient H: Confidence Intervals 
        and the Advantages and Drawbacks of Interpretive Guidelines. 
        \emph{Measurement and Evaluation in Counseling and Development, 57(2),} 93105.}
\cr\cr {McNeish, D. (2018). Thanks coefficient alpha, we'll take it from here. 
	      \emph{Psychological Methods, 23(3),} 412433.}
\cr\cr {Revelle, W. (2026). psych: Procedures for Psychological, Psychometric, 
       and Personality Research. R package version 2.6.5, 
       https://CRAN.R-project.org/package=psych}
\cr\cr {Revelle, W., & Condon, D. M. (2019). Reliability from alpha to omega: A tutorial. 
	      \emph{Psychological Assessment, 31(12),} 13951411.}
\cr\cr {Waller, N. G. (2018) Direct Schmid-Leiman Transformations and Rank-Deficient 
        Loadings Matrices. \emph{Psychometrika, 83(4),} 858870.}
}
\author{Brian P. O'Connor }
\examples{
OMEGA(data_RSE, Nfactors = 3, display = 1)
\donttest{
}
}