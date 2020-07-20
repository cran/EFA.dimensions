\name{data_RSE}
\alias{data_RSE}
\docType{data}
\title{Item-level dataset for the Rosenberg Self-Esteem scale}
\description{A data frame with 300 observations on the 10 items from the Rosenberg Self-Esteem scale.}
\usage{data(data_RSE)}
\examples{
head(data_RSE)

\donttest{
# MAP test on the Rosenberg Self-Esteem Scale (RSE) data
MAP(data_RSE, corkind='pearson', verbose=TRUE)

# parallel analysis of the Rosenberg Self-Esteem Scale (RSE) data
RAWPAR(data_RSE, extract='PCA', Ndatasets=100, percentile=95,
       corkind='pearson', verbose=TRUE)
}
}
