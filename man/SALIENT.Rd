\name{SALIENT}
\alias{SALIENT}
\title{Salient loadings criterion for determining the number of factors.}
\description{Salient loadings criterion for determining the number of factors,
	as recommended by Gorsuch. Factors are retained when they consist of a specified minimum 
	number (or more) variables that have a specified minimum (or higher) loading value.}
\usage{SALIENT(data, salvalue, numsals, corkind, verbose)}
\arguments{
 \item{data}{
 \code{}An all-numeric dataframe where the rows are cases & the columns are the variables,
        or a correlation matrix with ones on the diagonal. The function internally determines 
        whether the data are a correlation matrix.}

 \item{salvalue}{
 \code{}The loading value that is considered salient. Default = .40}

 \item{numsals}{
 \code{}The required number of salient loadings for a factor. Default = 3}

 \item{corkind}{
 \code{}The kind of correlation matrix to be used if data is not a correlation matrix. 
        The options are 'pearson', 'kendall', 'spearman', and 'polychoric'. Required only if
        the entered data is not a correlation matrix.}
                  
  \item{verbose}{
  \code{}Should detailed results be displayed in console? \cr The options are: TRUE (default) or FALSE.}
}
\value{The number of factors according to the salient loadings criterion.}
\references{
	      {Gorsuch, R. L. (1997a). Exploratory factor analysis: Its role in item analysis. 
	       \emph{Journal of Personality Assessment, 68,} 532-560.}
  \cr\cr {Boyd, K. C. (2011). Factor analysis. In M. Stausberg & S. Engler (Eds.), \emph{The Routledge
   	       Handbook of Research Methods in the Study of Religion} (pp. 204-216). New York: Routledge.}
}
\author{Brian P. O'Connor }
\examples{
# test on the Harman correlation matrix
SALIENT(data_Harman, salvalue=.4, numsals=3, corkind='pearson', verbose=TRUE)

\donttest{
# test on the Rosenberg Self-Esteem Scale (RSE) using Pearson correlations
SALIENT(data_RSE, salvalue=.4, numsals=3, corkind='pearson', verbose=TRUE)

# test on the Rosenberg Self-Esteem Scale (RSE) using polychoric correlations
SALIENT(data_RSE, salvalue=.4, numsals=3, corkind='polychoric', verbose=TRUE)
}
}
