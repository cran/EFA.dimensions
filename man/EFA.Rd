\name{EFA}
\alias{EFA}
\title{Exploratory factor analysis}
\description{Exploratory factor analysis with multiple options for factor extraction and rotation}
\usage{EFA(data, extraction = 'paf', corkind='pearson', Nfactors=NULL, Ncases=NULL, iterpaf=100, 
    rotation='promax', ppower = 3, verbose=TRUE)}
\arguments{
  \item{data}{
  \code{}An all-numeric dataframe where the rows are cases & the columns are 
  the variables, or a correlation matrix with ones on the diagonal.The function 
  internally determines whether the data are a correlation matrix.}

  \item{extraction}{
  \code{}The factor extraction method for the analysis. 
  The options are 'paf' (the default), 'ml', 'image', 'minres', 'uls', 'ols', 'wls', 'gls', 
  'alpha', and 'fullinfo'. }

  \item{corkind}{
  \code{}The kind of correlation matrix to be used if data is not a correlation matrix. 
  The options are 'pearson', 'kendall', 'spearman', 'gamma', and 'polychoric'. Required 
  only if the entered data is not a correlation matrix.}

  \item{Nfactors}{
  \code{}The number of factors to extract. If not specified, then the EMPKC procedure will be used
  to determine the number of factors.}

   \item{Ncases}{
  \code{}The number of cases. Required only if data is a correlation matrix.}

  \item{iterpaf}{
  \code{}The maximum number of iterations for paf.}

 \item{rotation}{
  \code{}The factor rotation method for the analysis. The orthogonal rotation options are: 
  'varimax' (the default), 'quartimax', 'bentlerT', 'equamax', 'geominT', 'bifactorT', 
  'entropy', and 'none'. The oblique rotation options are: 'promax' (the default), 'quartimin', 
  'oblimin', 'oblimax', 'simplimax', 'bentlerQ', 'geominQ', 'bifactorQ', and 'none'.}

  \item{ppower}{ 
  \code{}The power value to be used in a promax rotation (required only if
         rotation = 'promax'). Suggested value: 3}
                  
  \item{verbose}{
  \code{}Should detailed results be displayed in console? TRUE (default) or FALSE}
}
\details{The factor extraction computations for the following methods are conducted
	using the psych package (Revelle, 2023): 'minres', 'uls', 'ols', 'wls', 'gls', and 'alpha'.
	
	The factor extraction computations for 'fullinfo' are conducted
	using the mirt package (Chalmers, 2012). Full-information methods are considered more
	appropriate for item-level data than other factor extraction methods (Wirth & Edwards, 2007).
		
	The factor rotation computations for the following methods are conducted
	using the GPArotation package (Bernaards & Jennrich, 2005, 2023): 'quartimax', 'bentlerT', 
	'geominT', 'bifactorT', 'entropy', 'quartimin', 'oblimin', 'oblimax', 'simplimax', 
	'bentlerQ', 'geominQ', and 'bifactorQ'.
	
	For factor extraction (see Mulaik, 2010, for a review):
	 	
	\itemize{		
    \item\strong{paf} is for principal axis factor analysis
	     
    \item\strong{ml} is for maximum likelihood factor analysis
	      
    \item\strong{image} is for image factor analysis
	   
    \item\strong{minres} is for a minimum residual factor analysis (Revelle, 2023)

    \item\strong{uls} is for an unweighted least squares factor analysis (Revelle, 2023)
	   
    \item\strong{ols} is for an ordinary least squares factor analysis (Revelle, 2023)
	   
    \item\strong{wls} is for a weighted least squares factor analysis (Revelle, 2023)
	   
    \item\strong{gls} is for a generalized weighted least squares factor analysis (Revelle, 2023)
	   
    \item\strong{alpha} is for an alpha factor analysis (Revelle, 2023)
	}
	
    For factor rotation (see Jennrich, 2018, for a review):
    
	\itemize{
		
    \item\strong{varimax} is an orthogonal rotation that maximizes the spread of loadings within factors,
    which facilitates the interpretation of factors
    
    \item\strong{quartimax} is an orthogonal rotation that maximizes the spread of loadings for each variable 
    across factors, which facilitates the interpretation of variables (Bernaards & Jennrich, 2023)
        
    \item\strong{bentlerT} is an orthogonal rotation based on Bentler's invariant pattern simplicity 
    criterion (Bernaards & Jennrich, 2023)

    \item\strong{equamax} is an orthogonal rotation from the Crawford-Ferguson family (Bernaards & Jennrich, 2023)
        
    \item\strong{geominT} is an orthogonal rotation (Bernaards & Jennrich, 2023)
        
    \item\strong{bifactorT} is an orthogonal Jennrich and Bentler bifactor rotation (Bernaards & Jennrich, 2023)
        
    \item\strong{entropy} is a minimum entropy orthogonal rotation (Bernaards & Jennrich, 2023)
        
    \item\strong{promax} is an oblique rotation
        
    \item\strong{quartimin} is an oblique rotation (Bernaards & Jennrich, 2023)
        
    \item\strong{oblimin} is an oblique rotation (Bernaards & Jennrich, 2023)
        
    \item\strong{oblimax} is an oblique rotation (Bernaards & Jennrich, 2023)
        
    \item\strong{simplimax} is an oblique rotation (Bernaards & Jennrich, 2023)
        
    \item\strong{bentlerQ} is an oblique rotation based on Bentler''s invariant pattern simplicity 
    criterion (Bernaards & Jennrich, 2023)
    
    \item\strong{geominQ} is an oblique rotation (Bernaards & Jennrich, 2023)
    
    \item\strong{bifactorQ} is an oblique Jennrich and Bentler biquartimin rotation (Bernaards & Jennrich, 2023)
  	}
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
	  {Bernaards, C. A., & Jennrich, R. I. (2005). Gradient Projection Algorithms and Software for Arbitrary 
	  	Rotation Criteria in Factor Analysis. \emph{Educational and Psychological Measurement, 65(5),} 676-696. 
	  	https://doi.org/10.1177/0013164404272507}
\cr\cr{Bernaards, C. A., & Jennrich, R. I. (2023). GPArotation: Gradient Projection Factor Rotation. 
	R package version 2023.3-1, https://CRAN.R-project.org/package=GPArotation}
\cr\cr {Chalmers, R. P. (2012). mirt: A Multidimensional Item Response Theory Package for the R Environment. 
   \emph{Journal of Statistical Software, 48(6),} 129. doi:10.18637/jss.v048.i06.}
\cr\cr{Jennrich, R. I. (2018). Rotation. In P. Irwing, T. Booth, & D. J. Hughes (Eds.), \emph{The Wiley handbook 
 	of psychometric testing: A multidisciplinary reference on survey, scale and test development} (pp. 279304). 
 	Wiley Blackwell. https://doi.org/10.1002/9781118489772.ch10}
\cr\cr{Mulaik, S. A. (2010). \emph{Foundations of factor analysis} (2nd ed.). Boca Raton, FL: Chapman 
	and Hall/CRC Press, Taylor & Francis Group.}
\cr\cr{Revelle, W. (2023). psych: Procedures for Psychological, Psychometric, and Personality Research. 
	R package version 2.3.6, https://CRAN.R-project.org/package=psych}
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
EFA(data=data_Harman, extraction = 'paf', Nfactors=2, Ncases=305, rotation='oblimin', verbose=TRUE)

\donttest{
# Rosenberg Self-Esteem scale items, using ml extraction & bifactorQ rotation
EFA(data=data_RSE, extraction = 'ml', corkind='polychoric', Nfactors=2, 
    rotation='bifactorQ', verbose=TRUE)

# Rosenberg Self-Esteem scale items, using full-information factor extraction
EFA(data=data_RSE, extraction = 'fullinfo', corkind='pearson', Nfactors=2, 
    rotation='none', verbose=TRUE)

# NEO-PI-R scales
EFA(data=data_NEOPIR, extraction = 'minres', corkind='pearson', Nfactors=5, 
    iterpaf=100, rotation='promax', ppower = 4, verbose=TRUE) 
}
}
