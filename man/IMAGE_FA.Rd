\name{IMAGE_FA}
\alias{IMAGE_FA}
\title{Image factor analysis}
\description{Image factor analysis}
\usage{IMAGE_FA(data, corkind, Nfactors, rotate, ppower, verbose)}
\arguments{
  \item{data}{
  \code{}An all-numeric dataframe where the rows are cases & the columns are the variables, 
  or a correlation matrix with ones on the diagonal.The function internally determines 
  whether the data are a correlation matrix.}

  \item{corkind}{
  \code{}The kind of correlation matrix to be used if data is not a correlation matrix. 
  The options are 'pearson', 'kendall', 'spearman', and 'polychoric'. Required only 
  if the entered data is not a correlation matrix.}

  \item{Nfactors}{
  \code{}The number of factors to extract.}

  \item{rotate}{
  \code{}The factor rotation method. The options are: 'promax', 'varimax', and 'none'.}

  \item{ppower}{ 
  \code{}The power value to be used in a promax rotation (required only if
         rotate = 'promax'). Suggested value: 3}
                  
  \item{verbose}{
  \code{}Should detailed results be displayed in console? \cr The options are: TRUE (default) or FALSE.}
}
\value{ A list with the following elements:
 \item{eigenvalues}{The eigenvalues}
 \item{loadingsNOROT}{The unrotated factor loadings}
 \item{loadingsROT}{The rotated factor loadings (for varimax rotation)}
 \item{structure}{The structure matrix (for promax rotation)}
 \item{pattern}{The pattern matrix (for promax rotation)}
 \item{correls}{The correlations between the factors (for promax rotation)}
}
\author{Brian P. O'Connor }
\examples{
IMAGE_FA(data_NEOPIR, corkind='pearson', Nfactors=5, rotate='varimax', ppower=3, verbose=TRUE)
}