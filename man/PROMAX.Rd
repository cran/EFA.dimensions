\name{PROMAX}
\alias{PROMAX}
\title{promax rotation}
\description{promax rotation}
\usage{PROMAX(loadings, ppower, verbose)}
\arguments{
  \item{loadings}{
  \code{}A loading matrix.}

  \item{ppower}{
  \code{}The exponent for the promax target matrix. 'ppower' must be 1 or greater. 
         '4' is a conventional value.}
                  
  \item{verbose}{
  \code{}Should detailed results be displayed in console? TRUE (default) or FALSE}
}
\details{This function uses the R built-in promax function and provides additional output.}
\value{A list with the following elements:
 \item{loadingsNOROT}{The unrotated loadings}
 \item{pattern}{The pattern matrix (for promax rotation)}
 \item{structure}{The structure matrix (for promax rotation)}
 \item{phi}{The correlations between the factors (for promax rotation)}
 \item{eigenvar}{Eigenvalues and factor proportions of variance (for promax rotation)}
 \item{cormat_reproduced}{The reproduced correlation matrix, based on the rotated loadings}
}
\author{Brian P. O'Connor}
\examples{
\donttest{
# the Harman (1967) correlation matrix
PCAoutput <- PCA(data_Harman, Nfactors = 2, Ncases=305, rotate='none', verbose=TRUE)
PROMAX(PCAoutput$loadingsNOROT, ppower = 4, verbose=TRUE)

# Rosenberg Self-Esteem scale items
PCAoutput <- PCA(data_RSE, corkind='polychoric', Nfactors = 2, rotate='none', verbose=TRUE)
PROMAX(PCAoutput$loadingsNOROT, ppower = 4, verbose=TRUE)

# NEO-PI-R scales
PCAoutput <- PCA(data_NEOPIR, corkind='pearson', Nfactors = 5, rotate='none', verbose=TRUE)
PROMAX(PCAoutput$loadingsNOROT, ppower = 4, verbose=TRUE)
}
}
