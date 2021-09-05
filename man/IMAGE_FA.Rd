\name{IMAGE_FA}
\alias{IMAGE_FA}
\title{Image factor analysis}
\description{Image factor analysis}
\usage{IMAGE_FA(data, corkind, Nfactors=NULL, Ncases=NULL, rotate, ppower, verbose)}
\arguments{
  \item{data}{
  \code{}An all-numeric dataframe where the rows are cases & the columns are the variables, 
  or a correlation matrix with ones on the diagonal.The function internally determines 
  whether the data are a correlation matrix.}

  \item{corkind}{
  \code{}The kind of correlation matrix to be used if data is not a correlation matrix. 
  The options are 'pearson', 'kendall', 'spearman', 'gamma', and 'polychoric'. Required only 
  if the entered data is not a correlation matrix.}

  \item{Nfactors}{
  \code{}The number of factors to extract.}

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
\details{
"Image analysis is a principal factor variant in its usual application. As
in the case of the principal axes with estimated communalities procedure,
it is a principal factor variant in the sense that, after appropriate
alterations are made to the correlation matrix, that matrix can be
submitted to a principal factor program to find the desired factors. The
program will then minimize the residuals of the particular matrix
submitted. The differences in image analysis and other principal factor
variants lie primarily in the alterations of the correlation matrix before
the factors are extracted. Image factors can also be extracted by maximum
likelihood procedures (Joreskog, 1969b)." (Gorsuch, 1974, p. 103)

Image analysis is a common factor analysis of the predictable variance of a
set of variables. Each variable, Xi, is regressed on all other variables to
obtain the predicted Xi values. The covariance matrix of the predicted Xis
is then factored to produce an image analysis. The image of a variable is
the portion of a variable which is predicted from other variables.

The present function is an implementation of the Harris (1962) variant of
this procedure (see also Velicer, 1974), which is a noniterative
approximation to canonical component analysis. The results are identical to
the image factor analysis results that are produced by SAS and SPSS.
}
\value{A list with the following elements:
 \item{totvarexplNOROT}{The eigenvalues and total variance explained}
 \item{totvarexplROT}{The rotation sums of squared loadings and total variance explained for the rotated loadings}
 \item{loadingsNOROT}{The unrotated factor loadings}
 \item{loadingsROT}{The rotated factor loadings (for varimax rotation)}
 \item{structure}{The structure matrix (for promax rotation)}
 \item{pattern}{The pattern matrix (for promax rotation)}
 \item{correls}{The correlations between the factors (for promax rotation)}
 \item{cormat_reproduced}{The reproduced correlation matrix, based on the rotated loadings}
}
\references{ 
       {Gorsuch, R. L. (1983). \emph{Factor analysis ( 2nd ed. ).} Hillsdale, N.J.: Erlbaum, 1983.}
\cr\cr {Guttman, L. (1953). Image theory for the structure of quantitative variates. 
	\emph{Psychometrika 18,} 277-296.}
\cr\cr {Harris, C. W. (1962). Some Rao-Guttman relationships. \emph{Psychometrika, 27,}  247-63.}
\cr\cr {Velicer, W. F. (1974). A comparison of the stability of factor analysis, principal 
	component analysis, and rescaled image analysis. \emph{Educational and Psychological 
		Measurement, 34(3),} 563-572.}
}
\author{Brian P. O'Connor }
\examples{
\donttest{
# the Harman (1967) correlation matrix
IMAGE_FA(data_Harman, Nfactors=2, Ncases=305, rotate='PROMAX', ppower = 4, verbose=TRUE)

# Rosenberg Self-Esteem scale items
IMAGE_FA(data_RSE, corkind='gamma', Nfactors=2, rotate='PROMAX', ppower = 4, verbose=TRUE)

# NEO-PI-R scales
IMAGE_FA(data_NEOPIR, corkind='pearson', Nfactors=5, rotate='PROMAX', ppower = 4, verbose=TRUE)
}
}
