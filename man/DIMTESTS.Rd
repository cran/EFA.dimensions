\name{DIMTESTS}
\alias{DIMTESTS}
\title{Tests for the number of factors}
\description{Conducts multiple tests for the number of factors
}
\usage{DIMTESTS(data, tests, corkind, Ncases, HULL_method, HULL_gof, HULL_cor_method,
          CD_cor_method, display)}
\arguments{
  \item{data}{
  \code{}An all-numeric dataframe where the rows are cases & the columns are the variables, 
         or a correlation matrix with ones on the diagonal. The function internally determines 
         whether the data are a correlation matrix.}

  \item{tests}{
  \code{}A vector of the names of the tests for the number of factors that should
         be conducted. The possibilities are CD, EMPKC, HULL, MAP, NEVALSGT1, 
         RAWPAR, SALIENT, SESCREE, SMT. If tests is not specified, then 
         tests = c('EMPKC', 'HULL', 'RAWPAR') is used as the default.}
         
  \item{corkind}{
  \code{}The kind of correlation matrix to be used if data is not a correlation matrix. 
         The options are 'pearson', 'kendall', 'spearman', 'gamma', and 'polychoric'. Required only if 
         the entered data is not a correlation matrix.}
                  
  \item{Ncases}{
  \code{}The number of cases. Required only if data is a correlation matrix.}

  \item{HULL_method}{
  \code{}From EFAtools: The estimation method to use. One of "PAF" (default), "ULS", or "ML", 
         for principal axis factoring, unweighted least squares, and maximum likelihood}

  \item{HULL_gof}{
  \code{}From EFAtools: The goodness of fit index to use. Either "CAF" (default), "CFI", or "RMSEA", 
         or any combination of them. If method = "PAF" is used, only the CAF can be used as 
         goodness of fit index. For details on the CAF, see Lorenzo-Seva, Timmerman, and Kiers (2011).}

  \item{HULL_cor_method}{
  \code{}From EFAtools: The kind of correlation matrix to be used for the Hull method analyses. 
         The options are 'pearson', 'kendall', and 'spearman'}

  \item{CD_cor_method}{
  \code{}From EFAtools: The kind of correlation matrix to be used for the CD method analyses. 
         The options are 'pearson', 'kendall', and 'spearman'}

  \item{display}{
  \code{}The results to be displayed in the console: 0 = nothing; 1 = only the # of 
         factors for each test; 2 (default) = detailed output for each test}
}
\details{This is a convenience function for tests for the number of factors.
	
	The HULL method option uses the HULL function (and its defaults) in the EFAtools package.

	From Auerswald & Moshagen (2019):

	"The Hull method (Lorenzo-Seva et al., 2011) is an approach based on the Hull 
	heuristic used in other areas of model selection (e.g., Ceulemans & Kiers, 2006). 
	Similar to nongraphical variants of Cattell's scree plot, the Hull method attempts 
	to find an elbow as justification for the number of common factors. However, instead 
	of using the eigenvalues relative to the number of factors, the Hull method relies 
	on goodness-of-fit indices relative to the model degrees of freedom of the proposed 
	model."
	

	The CD (comparison data) method option uses the CD function (and its defaults) in the 
	EFAtools package. The CD method can only be conducted on raw data and not on 
	correlation matrices.

	From Auerswald & Moshagen (2019):

	"Ruscio and Roche (2012) suggested an approach that finds the number of factors 
	by determining the solution that reproduces the pattern of eigenvalues best 
	(comparison data, CD). CD takes previous factors into account by generating 
	comparison data of a known factorial structure in an iterative procedure. 
	Initially, CD compares whether the simulated comparison data with one underlying 
	factor (j = 1) reproduce the pattern of empirical eigenvalues significantly 
	worse compared with a two-factor solution (j + 1). If this is the case, CD 
	increases j until further improvements are nonsignificant or a preset maximum 
	of factors is reached."
	
	"No single extraction criterion performed best for every factor model. In 
	unidimensional and orthogonal models, traditional PA, EKC, and Hull 
	consistently displayed high hit rates even in small samples. Models with 
	correlated factors were more challenging, where CD and SMT outperformed 
	other methods, especially for shorter scales. Whereas the presence of 
	cross-loadings generally increased accuracy, non-normality had virtually 
	no effect on most criteria. We suggest researchers use a combination of 
	SMT and either Hull, the EKC, or traditional PA, because the number of 
	factors was almost always correctly retrieved if those methods converged. 
	When the results of this combination rule are inconclusive, traditional 
	PA, CD, and the EKC performed comparatively well. However, disagreement 
	also suggests that factors will be harder to detect, increasing sample 
	size requirements to N >= 500."
		

	The recommended tests for the number of factors are: EMPKC, HULL, and RAWPAR.
	The MAP test is also recommended for principal components analyses.
	Other possible methods (e.g., NEVALSGT1, SALIENT, SESCREE) are less 
	well-validated and are included for research purposes.
}
\value{A list with the following elements:
 \item{dimtests}{A matrix with the DIMTESTS results}
 
 \item{NfactorsDIMTESTS}{The number of factors according to the first test method 
 	                     specified in the "tests" vector}
}
\references{
{Auerswald, M., & Moshagen, M. (2019). How to determine the number of factors to 
	retain in exploratory factor analysis: A comparison of extraction methods under 
	realistic conditions. \emph{Psychological Methods, 24(4),} 468-491.}
\cr\cr {Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. (2011). The Hull method 
	for selecting the number of common factors. \emph{Multivariate Behavioral 
	Research, 46(2),} 340-364.}
\cr\cr {O'Connor, B. P. (2000). SPSS and SAS programs for determining 
  the number of components using parallel analysis and Velicer's 
  MAP test. \emph{Behavior Research Methods, Instrumentation, and
  Computers, 32,} 396-402.}
\cr\cr {Ruscio, J., & Roche, B. (2012). Determining the number of factors to retain 
	in an exploratory factor analysis using comparison data of known factorial 
	structure. \emph{Psychological Assessment, 24,} 282292. doi: 10.1037/a0025697}
\cr\cr {Zwick, W. R., & Velicer, W. F. (1986). Comparison of five rules for determining  
	the number of components to retain. \emph{Psychological Bulletin, 99,} 432-442.}
}
\author{Brian P. O'Connor }
\examples{
\donttest{
# the Harman (1967) correlation matrix
DIMTESTS(data_Harman, tests = c('EMPKC','HULL','RAWPAR'), corkind='pearson', 
                                Ncases = 305, display=2)

# Rosenberg Self-Esteem scale items, all possible DIMTESTS
DIMTESTS(data_RSE, 
         tests = c('CD','EMPKC','HULL','MAP','NEVALSGT1','RAWPAR','SALIENT','SESCREE','SMT'), 
      corkind='pearson', display=2)
  
# Rosenberg Self-Esteem scale items, using polychoric correlations
DIMTESTS(data_RSE, corkind='polychoric', display=2)

# NEO-PI-R scales
DIMTESTS(data_NEOPIR, tests = c('EMPKC','HULL','RAWPAR','NEVALSGT1'), display=2)
}
}
