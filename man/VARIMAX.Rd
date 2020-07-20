\name{VARIMAX}
\alias{VARIMAX}
\title{varimax rotation}
\description{varimax rotation}
\usage{VARIMAX(loadings, verbose)}
\arguments{
  \item{loadings}{
  \code{}A loading matrix.}
                  
  \item{verbose}{
  \code{}Should detailed results be displayed in console? \cr The options are: TRUE (default) or FALSE.}
}
\value{The varimax-rotated loadings}
\author{Brian P. O'Connor }
\examples{
\donttest{
loadings <- PCA(data_NEOPIR, corkind='pearson', Nfactors = 5, rotate='none', verbose=TRUE)

VARIMAX(loadings, verbose=TRUE)
}
}
