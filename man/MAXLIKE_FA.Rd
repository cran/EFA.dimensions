\name{MAXLIKE_FA}
\alias{MAXLIKE_FA}
\title{Maximum likelihood factor analysis}
\description{Maximum likelihood factor analysis}
\usage{MAXLIKE_FA(data, corkind, Nfactors, tolerml, iterml, rotate, ppower, verbose)}
\arguments{
  \item{data}{
  \code{}An all-numeric dataframe where the rows are cases & the columns are the 
  variables, or a correlation matrix with ones on the diagonal.The function internally 
  determines whether the data are a correlation matrix.}

  \item{corkind}{
  \code{}The kind of correlation matrix to be used if data is not a correlation matrix. 
  The options are 'pearson', 'kendall', 'spearman', and 'polychoric'. Required only 
  if the entered data is not a correlation matrix.}

  \item{Nfactors}{
  \code{}The number of factors to extract.}

  \item{tolerml}{
  \code{}The tolerance level.}

  \item{iterml}{
  \code{}The maximum number of iterations.}

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
MAXLIKE_FA(data_RSE, corkind='pearson', Nfactors = 2, 
           tolerml = .001, iterml = 50, rotate='promax', ppower=3, verbose=TRUE)
}
