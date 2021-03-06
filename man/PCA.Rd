\name{PCA}
\alias{PCA}
\title{Principal components analysis}
\description{Principal components analysis}
\usage{PCA(data, corkind, Nfactors=NULL, Ncases=NULL, rotate, ppower, verbose)}
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
  \code{}The number of components to extract.}

   \item{Ncases}{
  \code{}The number of cases. Required only if data is a correlation matrix.}

 \item{rotate}{
  \code{}The factor rotation method. The options are: 'PROMAX', , and 'none'.}

  \item{ppower}{ 
  \code{}The power value to be used in a promax rotation (required only if
         rotate = 'PROMAX'). Suggested value: 3}
                  
  \item{verbose}{
  \code{}Should detailed results be displayed in console? TRUE (default) or FALSE}
}
\value{A list with the following elements:
 \item{eigenvalues}{The eigenvalues}
 \item{loadingsNOROT}{The unrotated factor loadings}
 \item{loadingsROT}{The rotated factor loadings (for varimax rotation)}
 \item{structure}{The structure matrix (for promax rotation)}
 \item{pattern}{The pattern matrix (for promax rotation)}
 \item{correls}{The correlations between the factors (for promax rotation)}
 \item{cormat_reproduced}{The reproduced correlation matrix, based on the rotated loadings}
 \item{fit_coefficients}{Model fit coefficients}
}
\author{Brian P. O'Connor }
\examples{
\donttest{
# the Harman (1967) correlation matrix
PCA(data_Harman, Nfactors=2, Ncases=305, rotate='PROMAX', ppower = 4, verbose=TRUE)

# Rosenberg Self-Esteem scale items
PCA(data_RSE, corkind='polychoric', Nfactors=2, rotate='PROMAX', ppower = 4, verbose=TRUE)

# NEO-PI-R scales
PCA(data_NEOPIR, corkind='pearson', Nfactors=5, rotate='PROMAX', ppower = 4, verbose=TRUE)
}
}
