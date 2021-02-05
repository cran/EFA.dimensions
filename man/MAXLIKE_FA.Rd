\name{MAXLIKE_FA}
\alias{MAXLIKE_FA}
\title{Maximum likelihood factor analysis}
\description{Maximum likelihood factor analysis}
\usage{MAXLIKE_FA(data, corkind, Nfactors=NULL, Ncases=NULL, rotate, ppower, verbose)}
\arguments{
  \item{data}{
  \code{}An all-numeric dataframe where the rows are cases & the columns are the 
  variables, or a correlation matrix with ones on the diagonal.The function internally 
  determines whether the data are a correlation matrix.}

  \item{corkind}{
  \code{}The kind of correlation matrix to be used if data is not a correlation matrix. 
  The options are 'pearson', 'kendall', 'spearman', 'gamma', and 'polychoric'. Required only 
  if the entered data is not a correlation matrix.}

  \item{Nfactors}{
  \code{}The number of factors to extract.}

  \item{Ncases}{
  \code{}The number of cases. Required only if data is a correlation matrix.}

  \item{rotate}{
  \code{}The factor rotation method. The options are: 'PROMAX', 'VARIMAX', and 'none'.}

  \item{ppower}{ 
  \code{}The power value to be used in a promax rotation (required only if
         rotate = 'PROMAX'). Suggested value: 3}
                  
  \item{verbose}{
  \code{}Should detailed results be displayed in console? TRUE (default) or FALSE}
}
\details{This function relies heavily on the R factanal function, and it uses
	the fa from the psych package when factanal produces an error.	
}
% The code for this function was adapted from Reyment, Joreskog, & Marcus (1996, p. 308)}
\value{A list with the following elements:
 \item{eigenvalues}{The eigenvalues}
 \item{loadingsNOROT}{The unrotated factor loadings}
 \item{loadingsROT}{The rotated factor loadings (for varimax rotation)}
 \item{structure}{The structure matrix (for promax rotation)}
 \item{pattern}{The pattern matrix (for promax rotation)}
 \item{correls}{The correlations between the factors (for promax rotation)}
 \item{cormat_reproduced}{The reproduced correlation matrix, based on the rotated loadings}
 \item{chisqMODEL}{The model chi square statistic}
 \item{dfMODEL}{The model degrees of freedom}
 \item{pvalue}{The model p-value}
 \item{fit_coefficients}{Model fit coefficients}
}
\references{ 
	  Reyment, R., Joreskog, K., & Marcus, L. F. (1996). \emph{Applied Factor Analysis 
       in the Natural Sciences}. Cambridge, MA: Cambridge University Press.
}
\author{Brian P. O'Connor }
\examples{
\donttest{
# the Harman (1967) correlation matrix
MAXLIKE_FA(data_Harman, Nfactors = 2, Ncases = 305,
           rotate='PROMAX', ppower = 4, verbose=TRUE)

# Rosenberg Self-Esteem scale items
MAXLIKE_FA(data_RSE, corkind='gamma', Nfactors = 2, 
           rotate='PROMAX', ppower = 4, verbose=TRUE)

# NEO-PI-R scales
MAXLIKE_FA(data_NEOPIR, corkind='pearson', Nfactors = 5, 
           rotate='PROMAX', ppower = 4, verbose=TRUE)
}
}


