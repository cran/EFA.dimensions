\name{PROMAX}
\alias{PROMAX}
\title{Promax rotation}
\description{Promax rotation}
\usage{PROMAX(loadings, ppower, verbose)}
\arguments{
  \item{loadings}{
  \code{}A loading matrix.}

  \item{ppower}{
  \code{}The exponent for the promax target matrix. 'ppower' must be 1 or greater. '4' is a conventional value.}
                  
  \item{verbose}{
  \code{}Should detailed results be displayed in console? \cr The options are: TRUE (default) or FALSE.}
}
\value{ A list with the following elements:
 \item{structure}{The structure matrix (for promax rotation)}
 \item{pattern}{The pattern matrix (for promax rotation)}
 \item{correls}{The correlations between the factors (for promax rotation)}
}
\author{Brian P. O'Connor}
\examples{
\donttest{
loadings <- PCA(data_NEOPIR, corkind='pearson', Nfactors = 5, rotate='none', verbose=TRUE)

PROMAX(loadings, ppower = 3, verbose=TRUE)
}
}
