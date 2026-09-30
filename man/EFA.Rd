\name{EFA}
\alias{EFA}
\title{Exploratory factor analysis}
\description{Exploratory factor analysis with multiple options for factor extraction and rotation}
\usage{EFA(data, Nfactors=NULL, extraction = 'paf', rotation='promax', corkind='pearson', 
    Ncases=NULL, iterpaf=100, ppower = 3, 
    GPA_options = list(delta = .01, 
                       epsilon = .00001, 
                       normalize = FALSE, 
                       maxit = 1000, 
                       randomStarts = 50),
    verbose=TRUE)}
\arguments{
  \item{data}{
  \code{}An all-numeric dataframe where the rows are cases & the columns are 
  the variables, or a correlation matrix with ones on the diagonal.The function 
  internally determines whether the data are a correlation matrix.}

  \item{Nfactors}{
  \code{}The number of factors to extract. If not specified, then the EMPKC procedure 
  will be used to determine the number of factors.}

  \item{extraction}{
  \code{}The factor extraction method for the analysis. 
  The options are 'paf' (the default), 'alpha', 'fullinfo', 'gls', 'image', 'minres',  
    'ml', 'ols', 'uls', and 'wls'. }

  \item{rotation}{
  \code{}The factor rotation method for the analysis.
  \cr\cr The orthogonal rotation options are: 
  'varimax' (the default), 'bentlerT', 'entropy', 
  'equamax', 'geominT', 'quartimax', and 'none'.
  \cr\cr The oblique rotation options are: 'promax' (the default), 'bentlerQ', 
  'geominQ', 'oblimin', 'oblimax', 'quartimin', 
  'simplimax', and 'none'.}

  \item{corkind}{
  \code{}The kind of correlation matrix to be used if data is not a correlation matrix. 
  The options are 'pearson', 'kendall', 'spearman', 'gamma', and 'polychoric'. Required 
  only if the entered data is not a correlation matrix.}

  \item{Ncases}{
  \code{}The number of cases. Required only if data is a correlation matrix.}

  \item{iterpaf}{
  \code{}The maximum number of iterations for paf.}

  \item{ppower}{ 
  \code{}The power value to be used in a promax rotation (required only if
         rotation = 'promax'). Suggested value: 3}
                  
  \item{GPA_options}{ 
  \code{}(optional) A list with options for when a GPArotation is used.
  See the Details section below.}
  
  \item{verbose}{
  \code{}Should detailed results be displayed in console? TRUE (default) or FALSE}
}
\details{The factor extraction computations for the following methods are conducted
	using the psych package (Revelle, 2026): 'alpha', 'gls', 'minres', 'ols', 'uls', and 'wls'.
	
	The factor extraction computations for 'fullinfo' are conducted
	using the mirt package (Chalmers, 2012). Full-information methods are considered more
	appropriate for item-level data than other factor extraction methods (Wirth & Edwards, 2007).
		
	The factor rotation computations for the following methods are conducted
	using the GPArotation package (Bernaards & Jennrich, 2005, 2026): 'bentlerQ', 
	'bentlerT', 'entropy', 'geominQ', 'geominT', 'oblimax', 'oblimin', 'quartimax', 
	'quartimin', and 'simplimax'.
	
	For factor extraction (see Mulaik, 2010, for a review):
	 	
	\itemize{		
    \item\strong{alpha} is for an alpha factor analysis

    \item\strong{gls} is for a generalized weighted least squares factor analysis
    
    \item\strong{image} is for image factor analysis
	   
    \item\strong{minres} is for a minimum residual factor analysis

    \item\strong{ml} is for maximum likelihood factor analysis
	      
    \item\strong{ols} is for an ordinary least squares factor analysis
	   
    \item\strong{paf} is for principal axis factor analysis
	     
    \item\strong{uls} is for an unweighted least squares factor analysis
	   
    \item\strong{wls} is for a weighted least squares factor analysis \cr
	}
	
    For factor rotation (see Jennrich, 2018, for a review):

	\itemize{
		
    \item\strong{bentlerQ} is an oblique rotation based on Bentler's invariant 
    pattern simplicity criterion. 
    
    \item\strong{bentlerT} is an orthogonal rotation based on Bentler's 
    invariant pattern simplicity criterion

    \item\strong{entropy} entropy is an orthogonal factor rotation 
    that pushes factor loadings to be either 
    strongly dominant (close to 1.0) or cleanly absent (close to 0.0), 
    reducing the overall informational "noise" of the matrix.
        
    \item\strong{equamax} is an orthogonal rotation designed as a 
    mathematical compromise between varimax and quartimax. 
        
    \item\strong{geominQ} is an oblique factor rotation method that 
    find a clean "simple structure" 
    even when the data contains highly complex variables (variables that 
    naturally load on multiple factors).
    
    \item\strong{geominT} is an orthogonal factor rotation that 
    combines the strict geometric constraint of uncorrelated factors with 
    Brownes flexible, product-based geomin complexity criterion. 
        
    \item\strong{oblimax} is an oblique factor rotation method that maximizes 
    the number of very high and very low (near-zero) factor loadings. 
        
    \item\strong{oblimin} is an oblique factor rotation that allows 
    factors to correlate freely to achieve the simplest possible loading structure. 
        
    \item\strong{promax} is a two-stage oblique factor rotation method. 
        
    \item\strong{quartimax} is an orthogonal factor rotation 
    designed to spread the explained variance evenly across the factors, 
    actively preventing a single general factor from dominating. 
        
    \item\strong{simplimax} is an oblique factor rotation that attempts to 
    recover complex or messy "simple structures" 
    where traditional continuous methods like oblimin or geomin fail.
        
    \item\strong{quartimin} is an oblique factor rotation that focuses 
    on simplifying the rows of the loading matrix.
        
    \item\strong{varimax} is an orthogonal 
    factor rotation that maximizes the variance of the squared factor loadings within 
    each individual factor column. \cr
  	}


  For \strong{GPA_options:}
  
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

See this web page for more detailed descriptions of the above 
extraction and rotation methods:

\url{https://oconnor-psych.ok.ubc.ca/vignettes/EFA_BIFACTOR_vignettes.html}
}
\value{A list with the following elements:
 \item{loadingsNOROT}{The unrotated factor loadings}
 \item{loadingsROT}{The rotated factor loadings}
 \item{pattern}{The pattern matrix}
 \item{structure}{The structure matrix}
 \item{phi}{The correlations between the factors}
 \item{varexplNOROT1}{The initial eigenvalues and total variance explained}
 \item{varexplNOROT2}{The eigenvalues and total variance explained after factor extraction (no rotation)}
 \item{varexplROT}{The rotation sums of squared loadings and total variance explained for the rotated loadings}
 \item{cormat_reprod}{The reproduced correlation matrix, based on the rotated loadings}
 \item{fit_coefs}{Model fit coefficients}
 \item{chisqMODEL}{The model chi squared}
 \item{dfMODEL}{The model degrees of freedom}
 \item{pvalue}{The model p-value}
 \item{chisqNULL}{The null model chi squared}
 \item{dfNULL}{The null model degrees of freedom}
 \item{communalities}{The unrotated factor solution communalities}
 \item{uniquenesses}{The unrotated factor solution uniquenesses}
}
\references{
	  {Bernaards, C. A., & Jennrich, R. I. (2005). Gradient Projection Algorithms  
	  	and Software for Arbitrary Rotation Criteria in Factor Analysis. 
	  	\emph{Educational and Psychological Measurement, 65(5),} 676-696. 
	  	https://doi.org/10.1177/0013164404272507}
\cr\cr{Bernaards, C. A., & Jennrich, R. I. (2026). GPArotation: Gradient Projection Factor Rotation. 
	R package version 2026.4-1, https://CRAN.R-project.org/package=GPArotation}
\cr\cr {Chalmers, R. P. (2012). mirt: A Multidimensional Item Response Theory Package for the R Environment. 
   \emph{Journal of Statistical Software, 48(6),} 129. doi:10.18637/jss.v048.i06.}
\cr\cr{Jennrich, R. I. (2018). Rotation. In P. Irwing, T. Booth, & D. J. Hughes (Eds.), \emph{The Wiley handbook 
 	of psychometric testing: A multidisciplinary reference on survey, scale and test development} (pp. 279304). 
 	Wiley Blackwell. https://doi.org/10.1002/9781118489772.ch10}
\cr\cr{Mulaik, S. A. (2010). \emph{Foundations of factor analysis} (2nd ed.). Boca Raton, FL: Chapman 
	and Hall/CRC Press, Taylor & Francis Group.}
\cr\cr{Revelle, W. (2026). psych: Procedures for Psychological, Psychometric, and Personality Research. 
	R package version 2.6.5, https://CRAN.R-project.org/package=psych}
\cr\cr {Sellbom, M., & Tellegen, A. (2019). Factor analysis in psychological assessment research: 
	Common pitfalls and recommendations. 
   \emph{Psychological Assessment, 31(12),} 1428-1441. https://doi.org/10.1037/pas0000623}
\cr\cr {Watts, A. L., Greene, A. L., Ringwald, W., Forbes, M. K., Brandes, C. M., Levin-Aspenson, 
	H. F., & Delawalla, C. (2023). Factor analysis in personality disorders research: Modern issues 
	and illustrations of practical recommendations. 
  \emph{Personality Disorders: Theory, Research, and Treatment, 14(1),} 105-117. 
  https://doi.org/10.1037/per0000581}
\cr\cr {Wirth, R. J., & Edwards, M. C. (2007). Item factor analysis: current approaches and future directions. 
  \emph{Psychological methods, 12(1),} 58-79. https://doi.org/10.1037/1082-989X.12.1.58}
}
\author{Brian P. O'Connor }
\examples{
# the Harman (1967) correlation matrix
EFA(data=data_Harman, extraction = 'paf', Nfactors=2, Ncases=305, rotation='oblimin')
\donttest{
# Rosenberg Self-Esteem scale items, using ml extraction & polychoric correlations
EFA(data=data_RSE, Nfactors=2, extraction = 'ml', corkind='polychoric')

# NEO-PI-R scales
EFA(data=data_NEOPIR, Nfactors=5, extraction='minres', rotation='promax') 
}
}
