\name{VARIMAX}
\alias{VARIMAX}
\title{varimax rotation}
\description{varimax rotation}
\usage{VARIMAX(loadings, normalize = TRUE, verbose)}
\arguments{
  \item{loadings}{
  \code{}A loading matrix.}
            
  \item{normalize}{
  \code{}Should Kaiser normalization be performed? If so the rows of x are re-scaled to 
  unit length before rotation, and scaled back afterwards. Default = TRUE.}
                 
  \item{verbose}{
  \code{}Should detailed results be displayed in console? TRUE (default) or FALSE}
}
\details{This function uses the R built-in varimax function and provides additional output.}
\value{A list with the following elements:
 \item{loadingsNOROT}{The unrotated loadings}
 \item{loadingsV}{The varimax-rotated loadings}
 \item{rotmatV}{The rotation matrix}
 \item{cormat_reproduced}{The reproduced correlation matrix, based on the rotated loadings}
}
\author{Brian P. O'Connor}
\examples{
\donttest{
# the Harman (1967) correlation matrix
PCAoutput <- PCA(data_Harman, Nfactors = 2, Ncases=305, rotate='none', verbose=TRUE)
VARIMAX(PCAoutput$loadingsNOROT, verbose=TRUE)

# Rosenberg Self-Esteem scale items
PCAoutput <- PCA(data_RSE, corkind='polychoric', Nfactors = 2, rotate='none', verbose=TRUE)
VARIMAX(PCAoutput$loadingsNOROT, verbose=TRUE)

# NEO-PI-R scales
PCAoutput <- PCA(data_NEOPIR, Nfactors = 5, rotate='none', verbose=TRUE)
VARIMAX(PCAoutput$loadingsNOROT, verbose=TRUE)
}
}
