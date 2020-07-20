\name{data_Harman}
\alias{data_Harman}
\docType{data}
\title{Correlation matrix from Harman (1967, p. 80).}
\description{The correlation matrix for eight physical variables for 305 cases from Harman (1967, p. 80).}
\usage{data(data_Harman)}
\references{Harman, H. H. (1967). \emph{Modern factor analysis (2nd. ed.).} Chicago: University of Chicago Press.}
\examples{
\donttest{
# MAP test on the Harman correlation matrix
MAP(data_Harman, verbose=TRUE)

# parallel analysis of the Harman correlation matrix
RAWPAR(data_Harman, extract='PCA', Ndatasets=100, percentile=95,
       Ncases=305, verbose=TRUE)
}
}
