\name{BIFACTOR}
\alias{BIFACTOR}
\title{Bifactor analysis}
\description{Bifactor analysis with multiple options for data entry and analysis}
\usage{BIFACTOR(loadings = NULL,
         rawdata = NULL,   
         cormat = NULL, Ncases = NULL, corkind = 'pearson', 
         Nfactors = 4,
         bifactor_kind = 'bifactorT', 
         EFA_options    = list(extraction = 'minres', 
                               rotation = 'oblimin'),
         LV_options     = list(group_keys = NULL,
                               estimator='MLR',
                               rotation = 'bigeomin',
                               resid_correls=NULL, 
                               LV_names=NULL,
                               ordered = FALSE),
         schmid_options = list(extraction = 'minres', 
                               rotation = 'oblimin'),
         GPA_options    = list(delta = .01, 
                               epsilon = .00001, 
                               normalize = FALSE, 
                               maxit = 1000, 
                               randomStarts = 50),
         min_loading = .2,
         verbose = TRUE)}
                 
\arguments{
  \item{loadings}{
  \code{}(optional) loadings can be either a bifactor loading matrix or an EFA
  loading matrix. See the Details and Examples sections below.}

  \item{rawdata}{
  \code{}(optional) An all-numeric dataframe where the rows are cases & the columns are 
  the variables. Required only when loadings and cormat are not provided.}

  \item{Nfactors}{
  \code{}(optional) The number of factors to extract, including the general factor.
  \cr Nfactors is = 1 (for the general factor) + the number of group factors.}

  \item{cormat}{
  \code{}(optional) A correlation matrix with ones on the diagonal.
  Required only when loadings and rawdata are not provided.}

  \item{Ncases}{
  \code{}(optional) The number of cases. Required only when cormat is provided.}

  \item{corkind}{
  \code{}(optional) The kind of correlation matrix to be used when rawdata is provided. 
  The options are 'pearson', 'kendall', 'spearman', 'gamma', and 'polychoric'.}

  \item{bifactor_kind}{
  \code{}(optional) Specify one of the following possible bifactor methods: 
  'bifactorT', 'bifactorQ', 'bigeominT', 'bigeominQ', 'SL', 'SLiD', 'DSL', 
  'CFA', 'ESEM', or 'none'.
  Use 'bifactor_kind = "none"' only when loadings is (already) a bifactor loading matrix.
  \cr Example: bifactor_kind = 'bifactorT'}

  \item{EFA_options}{
  \code{}(optional) A list with EFA options when rawdata is provided. The list elements
  include values for 'extraction' and 'rotation'.
  \cr\cr The possibilities for extraction are: 'minres'(the default), 'alpha', 'fullinfo', 
  'gls', 'image', 'ml', 'ols', 'paf', 'uls', and 'wls'.
  \cr\cr The possibilities for rotation are:  
  'varimax', 'bentlerT', 'entropy', 'equamax', 'geominT', 'quartimax', 
  'promax', 'bentlerQ', 'geominQ', 'oblimin' (the default, 'oblimax', 'quartimin', 
  'simplimax', and 'none'.}
 
  \item{LV_options}{
  \code{}(optional) A list with cfa options when bifactor_kind = 'CFA'. 
  The possible list elements are 'group_keys', 'N_group_factors',
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
  \code{}The minimum value of a group factor loading for an item to be
  considered to have a non-negligible contribution to a group factor.
  min_loading only plays a role in the computations for PUC and for
  some group factor statistics.}
                  
  \item{verbose}{
  \code{}Should detailed results be displayed in console? TRUE (default) or FALSE}
}
\details{
  The bifactor analyses can be conducted for three possible data input scenarios.

	\itemize{
    \item
  When "loadings" is provided and when it is (already) a bifactor loading matrix, 
  then bifactor rotations are not conducted and the function merely provides a variety
  of bifactor model, factor, and item statistics. "bifactor_kind" must be set
  to "none" for the function to recognize that loadings is (already) a 
  bifactor loading matrix.

    \item
  When "loadings" is provided and when it consists of EFA loadings (not bifactor 
  loadings), then the function will conduct the requested kind of bifactor
  analysis and provide the bifactor loadings and statistics.

    \item
  if either "rawdata" or "cormat" are provided, the function will conduct the 
  preliminary EFA, the bifactor analyses, and produce the bifactor loadings and  
  statistics. "loadings" must be NULL for this to occur.
  }

  The factor extraction computations for the following methods are conducted
	using the psych package (Revelle, 2026): 'alpha', 'gls', 'minres', 'ols', 'uls', and 'wls'.
	
	The factor extraction computations for 'fullinfo' are conducted
	using the mirt package (Chalmers, 2012). Full-information methods are considered more
	appropriate for item-level data than other factor extraction methods (Wirth & Edwards, 2007).
		
	The factor rotation computations for the following methods are conducted
	using the GPArotation package (Bernaards & Jennrich, 2005, 2026): 'bentlerQ', 
	'bentlerT', 'entropy', 'geominQ', 'geominT', 'oblimax', 'oblimin', 'quartimax', 
	'quartimin', and 'simplimax'.
	
	The 'bifactorT', 'bifactorQ', 'bigeominQ', and 'bigeominQ' computations 
	are conducted using the GPArotation package (Bernaards & Jennrich, 2005, 2026). 
	The 'SL' and 'DSL' bifactor computations 
	are conducted using the psych package (Revelle, 2026). The 'SLiD' bifactor  
	computations are conducted using code from Garcia-Garzon et al. (2021).	
	
	The only possible extraction methods for 'SL', 'SLiD', and 'DSL' bifactor 
	computations are 'paf', 'minres', and 'ml'.	\cr
	
  For reviews of bifactor analyses, see Bornovalova (2020), 
  Garcia-Garzon et al. (2021), & Reise et al. (2018).

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

  When bifactor_kind = 'CFA' or 'ESEM', then LV_options is a list with information for the
  theseanalyses. The possible list elements are 'group_keys', 
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

See this web page for more detailed descriptions of the above extraction, 
rotation, and bifactor methods:

\url{https://oconnor-psych.ok.ubc.ca/vignettes/EFA_BIFACTOR_vignettes.html}
}
\value{A list with the bifactor loadings and statistics, including omegas, ECV,
ARPB, FD, coef_H, and rmsr.
}
\references{
	  {Bernaards, C. A., & Jennrich, R. I. (2005). Gradient Projection Algorithms  
	  	and Software for Arbitrary Rotation Criteria in Factor Analysis. 
	  	\emph{Educational and Psychological Measurement, 65(5),} 676-696. 
	  	https://doi.org/10.1177/0013164404272507}
\cr\cr{Bernaards, C. A., & Jennrich, R. I. (2026). GPArotation: Gradient Projection  
	Factor Rotation. R package version 2026.4-1, 
	https://CRAN.R-project.org/package=GPArotation}
\cr\cr {Bornovalova, M. A., Choate, A. M., Fatimah, H., Petersen, K. J., & 
  Wiernik, B. M. (2020). Appropriate Use of Bifactor Analysis in Psychopathology 
  Research: Appreciating Benefits and Limitations. 
  \emph{Biological psychiatry, 88(1),} 1827.}
\cr\cr {Chalmers, R. P. (2012). mirt: A Multidimensional Item Response Theory 
  Package for the R Environment. 
   \emph{Journal of Statistical Software, 48(6),} 129. doi:10.18637/jss.v048.i06.}
\cr\cr{Garcia-Garzon, E., Abad, F. J., & Garrido, L. E. (2021). On omega  
        hierarchical estimation: A Comparison of Exploratory Bi-Factor Analysis
        Algorithms. \emph{Multivariate Behavioral Research, 56(1),} 101-119.}
\cr\cr{Giordano, C., Ones, D. S., Waller, N. G., & Stanek, K. C. (2020). 
  Exploratory bifactor measurement models in vocational behavior research. 
  \emph{Journal of Vocational Behavior, 120,} Article 103430.}
\cr\cr{Jennrich, R. I. (2018). Rotation. In P. Irwing, T. Booth, & D. J. Hughes (Eds.), 
  \emph{The Wiley handbook of psychometric testing: A multidisciplinary reference 
  on survey, scale and test development} (pp. 279304). 
 	Wiley Blackwell. https://doi.org/10.1002/9781118489772.ch10}
\cr\cr{Markon K. E. (2019). Bifactor and Hierarchical Models: Specification, 
  Inference, and Interpretation. \emph{Annual review of clinical psychology, 15,} 5169.}
\cr\cr{Mulaik, S. A. (2010). \emph{Foundations of factor analysis} (2nd ed.). Boca Raton, FL: Chapman 
	and Hall/CRC Press, Taylor & Francis Group.}
\cr\cr {Reise, S.P., Bonifay, W., & Haviland, M.G. (2018). Bifactor Modelling and 
  the Evaluation of Scale Scores. In P. Irwing, T. Booth& D.J. Hughes (eds.),
  \emph{The Wiley Handbook of Psychometric Testing}. John Wiley & Sons.}
\cr\cr{Revelle, W. (2026). psych: Procedures for Psychological, Psychometric, and 
	Personality Research. R package version 2.6.5, https://CRAN.R-project.org/package=psych}
\cr\cr {Sellbom, M., & Tellegen, A. (2019). Factor analysis in psychological assessment research: 
	Common pitfalls and recommendations. 
   \emph{Psychological Assessment, 31(12),} 1428-1441. https://doi.org/10.1037/pas0000623}
\cr\cr {Watts, A. L., Greene, A. L., Ringwald, W., Forbes, M. K., Brandes, C. M.,  
	Levin-Aspenson, H. F., & Delawalla, C. (2023). Factor analysis in personality 
	disorders research: Modern issues and illustrations of practical recommendations. 
  \emph{Personality Disorders: Theory, Research, and Treatment, 14(1),} 105-117. 
  https://doi.org/10.1037/per0000581}
\cr\cr {Wirth, R. J., & Edwards, M. C. (2007). Item factor analysis: current 
  approaches and future directions. \emph{Psychological methods, 12(1),} 58-79. 
  https://doi.org/10.1037/1082-989X.12.1.58}
}
\author{Brian P. O'Connor }
\examples{
# using only EFA loadings as input
efa_loadings <- EFA(data=data_RSE, Nfactors=3, extraction = 'ml', verbose=FALSE)$pattern

BIFACTOR(loadings = efa_loadings)  

\donttest{
# when a correlation matrix is provided, but not loadings or rawdata
BIFACTOR(cormat = cor(data_RSE), Ncases = 800, Nfactors = 4,
         bifactor_kind = 'SL', 
         schmid_options = list(extraction = 'minres', 
                               rotation = 'oblimin')) 
  
# when rawdata is provided, but not loadings or cormat
BIFACTOR(rawdata = data_RSE, Nfactors = 3,
         bifactor_kind = 'bifactorT', 
         EFA_options = list(extraction = 'minres', rotation = 'oblimin')) 
         
# confirmatory bifactor analyses
BIFACTOR(rawdata = data_RSE, 
         bifactor_kind = 'CFA', 
         LV_options = list(group_keys = c(1, 1, 2, 1, 2, 1, 1, 2, 2, 2), 
                           N_group_factors = NULL,
                           estimator = 'WLSMV',
                           rotation = 'bigeomin',
                           resid_correls = NULL,
                           LV_names=NULL,
                           ordered = TRUE),
         min_loading = .0001) 
}
# See this web page for additional Examples:

# https://oconnor-psych.ok.ubc.ca/vignettes/EXAMPLES_vignettes.html
}
