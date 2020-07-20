\name{SCREE_PLOT}
\alias{SCREE_PLOT}
\title{Scree plot of eigenvalues}
\description{Produces a scree plot of eigenvalues for raw data or for a correlation matrix.}
\usage{SCREE_PLOT(data, corkind)}
\arguments{
  \item{data}{
  \code{}An all-numeric dataframe where the rows are cases & the columns are 
  the variables, or a correlation matrix with ones on the diagonal.The function 
  internally determines whether the data are a correlation matrix.}

  \item{corkind}{
  \code{}The kind of correlation matrix to be used if data is not a correlation matrix. 
  The options are 'pearson', 'kendall', 'spearman', and 'polychoric'. Required 
  only if the entered data is not a correlation matrix.}

}
\author{Brian P. O'Connor }
\examples{
\donttest{
SCREE_PLOT(data_Field, corkind='pearson')

SCREE_PLOT(data_RSE, corkind='polychoric')
}
}
