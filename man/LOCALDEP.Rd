\name{LOCALDEP}
\alias{LOCALDEP}
\title{Local independence}
\description{Provides the residual correlations after partialling the first component 
	out of a correlation matrix. Item response theory models are based on the assumption that the items 
	display local independence. The latent trait is presumed to be responsible for 
	the associations between the items. Once the latent trait is partialled out, 
	the residual correlations between pairs of items should be negligible. Local 
	dependence exists when there is additional systematic covariance among the items. 
	It can occur when pairs of items have highly similar content or between sequentially 
	presented items in a test. Local dependence distorts IRT parameter estimates, 
	it can artificially increase scale information, and it distorts the latent trait, 
	which becomes too heavily defined by the locally dependent items. The LOCALDEP 
	function partials out the first component (not the IRT latent trait) from a 
	correlation matrix. Examining the residual correlations is a preliminary, 
	exploratory method of determining whether local dependence exists. The function 
	also displays the number of residual correlations that are >= a range of values.}
\usage{LOCALDEP(data, corkind, verbose)}
\arguments{
  \item{data}{
  \code{}An all-numeric dataframe where the rows are cases & the columns are the variables, 
  or a correlation matrix with ones on the diagonal. The function internally determines 
  whether the data are a correlation matrix.}

  \item{corkind}{
  \code{}The kind of correlation matrix to be used if data is not a correlation matrix. 
  The options are 'pearson', 'kendall', 'spearman', and 'polychoric'. Required only if 
  the entered data is not a correlation matrix.}
                  
  \item{verbose}{
  \code{}Should detailed results be displayed in console? \cr The options are: TRUE (default) or FALSE.}
}
\value{ A list with the following elements:
 \item{correlations}{The correlation matrix}
 \item{residcor}{The residualized correlation matrix}
}
\author{Brian P. O'Connor}
\examples{
# Residual correlations for the Rosenberg Self-Esteem Scale (RSE)
LOCALDEP(data_RSE, corkind = 'pearson', verbose=TRUE)
}
