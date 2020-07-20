\name{CONGRUENCE}
\alias{CONGRUENCE}
\title{Factor solution congruence}
\description{Aligns two factor loading matrices and computes the 
	factor solution congruence and the root mean square residual.}
\usage{CONGRUENCE(target, loadings, verbose)}
\arguments{
  \item{target}{
  \code{}The target loading matrix.}

  \item{loadings}{
  \code{}The loading matrix that will be aligned with the target.}
                  
  \item{verbose}{
  \code{}Should detailed results be displayed in console? \cr The options are: TRUE (default) or FALSE.}
}
\details{
The function first searches for the alignment of the factors from the two loading 
matrices that has the highest factor solution congruence. It then aligns the factors 
in "loadings" with the factors in "target" without changing the loadings. The alignment 
is based solely on the positions and directions of the factors. The function then 
produces the Tucker-Wrigley-Neuhaus factor solution congruence coefficient as an 
index of the degree of similarity between between the aligned loading matrices 
(see Guadagnoli & Velicer, 1991; and ten Berge, 1986, for reviews).
}
\value{ A list with the following elements:
 \item{rcBefore}{The factor solution congruence before factor alignment}
 \item{rcAfter}{The factor solution congruence after factor alignment}
 \item{rcFactors}{The congruence for each factor}
 \item{rmsr}{The root mean square residual}
 \item{residmat}{The residual matrix}
 \item{loadingsNew}{The aligned loading matrix}
}
\references{
{Guadagnoli, E., & Velicer, W. (1991). A comparison of pattern matching indices. 
	\emph{Multivariate Behavior Research, 26,} 323-343.}
\cr\cr{ten Berge, J. M. F. (1986). Some relationships between descriptive comparisons 
	of components from different studies. \emph{Multivariate Behavioral Research, 21,} 29-40.}
}
\author{Brian P. O'Connor}
\examples{
# RSE data
loadings <- PCA(data_RSE[1:150,],   corkind='pearson', Nfactors = 3, 
                rotate='varimax', verbose=FALSE)
target   <- PCA(data_RSE[151:300,], corkind='pearson', Nfactors = 3, 
                rotate='varimax', verbose=FALSE)
CONGRUENCE(target$loadingsROT, loadings$loadingsROT, verbose=TRUE)

\donttest{
# NEO-PI-R data
loadings <- PCA(data_NEOPIR[1:500,], corkind='pearson', Nfactors = 3, 
                rotate='varimax', verbose=FALSE)
target <- PCA(data_NEOPIR[501:1000,], corkind='pearson', Nfactors = 3, 
              rotate='varimax', verbose=FALSE)
CONGRUENCE(target$loadingsROT, loadings$loadingsROT, verbose=TRUE)
}
}
