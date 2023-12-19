\name{PCA}
\alias{PCA}
\title{Principal components analysis}
\description{Principal components analysis}
\usage{PCA(data, corkind='pearson', Nfactors=NULL, Ncases=NULL, rotation='promax', 
	 ppower=3, verbose=TRUE, rotate)}
\arguments{
  \item{data}{
  \code{}An all-numeric dataframe where the rows are cases & the columns are 
  the variables, or a correlation matrix with ones on the diagonal.The function 
  internally determines whether the data are a correlation matrix.}

  \item{corkind}{
  \code{}The kind of correlation matrix to be used if data is not a correlation matrix. 
  The options are 'pearson', 'kendall', 'spearman', 'gamma', and 'polychoric'. Required 
  only if the entered data is not a correlation matrix.}

  \item{Nfactors}{
  \code{}The number of components to extraction. If not specified, then the EMPKC procedure will be used
  to determine the number of components.}

   \item{Ncases}{
  \code{}The number of cases. Required only if data is a correlation matrix.}

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
  
  \item{rotate}{ 
  \code{}(Deprecated.) Use 'rotation' instead.}
}
\value{A list with the following elements:
 \item{loadingsNOROT}{The unrotated factor loadings}
 \item{loadingsROT}{The rotated factor loadings}
 \item{pattern}{The pattern matrix}
 \item{structure}{The structure matrix}
 \item{phi}{The correlations between the factors}
 \item{varexplNOROT1}{The initial eigenvalues and total variance explained}
 \item{varexplROT}{The rotation sums of squared loadings and total variance explained for the rotated loadings}
 \item{cormat_reprod}{The reproduced correlation matrix, based on the rotated loadings}
 \item{fit_coeffs}{Model fit coefficients}
 \item{communalities}{The unrotated factor solution communalities}
 \item{uniquenesses}{The unrotated factor solution uniquenesses}
}
\author{Brian P. O'Connor }
\examples{
\donttest{
# the Harman (1967) correlation matrix
PCA(data_Harman, Nfactors=2, Ncases=305, rotation='oblimin', verbose=TRUE)

# Rosenberg Self-Esteem scale items
PCA(data_RSE, corkind='polychoric', Nfactors=2, rotation='bifactorQ', verbose=TRUE)

# NEO-PI-R scales
PCA(data_NEOPIR, corkind='pearson', Nfactors=5, rotation='promax', ppower = 4, verbose=TRUE)
}
}
