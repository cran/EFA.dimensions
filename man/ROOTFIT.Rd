\name{ROOTFIT}
\alias{ROOTFIT}
\title{Factor fit coefficients}
\description{A variety of fit coefficients for the possible N-factor solutions 
	in exploratory factor analysis}
\usage{ROOTFIT(data, corkind, Ncases, extract, verbose)}
\arguments{
  \item{data}{
  \code{}An all-numeric dataframe where the rows are cases & the columns are the 
  variables, or a correlation matrix with ones on the diagonal.The function 
  internally determines whether the data are a correlation matrix.}

  \item{corkind}{
  \code{}The kind of correlation matrix to be used if data is not a correlation 
  matrix. The options are 'pearson', 'kendall', 'spearman', and 'polychoric'. Required 
  only if the entered data is not a correlation matrix.}

  \item{Ncases}{
  \code{}The number of cases upon which a correlation matrix is based. Required only 
  if data is a correlation matrix.}

  \item{extract}{ 
  \code{}The factor extraction method. The options are:
  	    'PAF' for principal axis / common factor analysis;
  	    'PCA' for principal components analysis.
  	    'ML' for maximum likelihood estimation.}
                  
  \item{verbose}{
  \code{}Should detailed results be displayed in console? \cr The options are: TRUE (default) or FALSE.}
}
\value{A list with eigenvalues & fit coefficients.}
\author{Brian P. O'Connor }
\examples{
# RSE data
ROOTFIT(data_RSE, corkind='pearson', extract='ML')
ROOTFIT(data_RSE, corkind='pearson', extract='PCA')

\donttest{
# NEO-PI-R data
ROOTFIT(data_NEOPIR, corkind='pearson', extract='ML')
ROOTFIT(data_NEOPIR, corkind='pearson', extract='PCA')
}
}
