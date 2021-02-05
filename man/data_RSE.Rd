\name{data_RSE}
\alias{data_RSE}
\docType{data}
\title{Item-level dataset for the Rosenberg Self-Esteem scale}
\description{A data frame with 300 observations on the 10 items from the Rosenberg Self-Esteem scale.}
\usage{data(data_RSE)}
\references{Rosenberg, M. (1965). \emph{Society and the adolescent self-image.} Princeton University Press.}
\examples{
\donttest{
# MAP test on the Rosenberg Self-Esteem Scale (RSE) item data
MAP(data_RSE, corkind='polychoric', verbose=TRUE)

# DIMTESTS on the Rosenberg Self-Esteem Scale (RSE) item data
DIMTESTS(data_RSE, tests = c('CD','EMPKC','HULL','RAWPAR','NEVALSGT1'), Ncases=1000, display=2)

# parallel analysis of the Rosenberg Self-Esteem Scale (RSE) item data
RAWPAR(data_RSE, factormodel='PCA', Ndatasets=100, percentile=95,
       corkind='pearson', verbose=TRUE)
}
}
