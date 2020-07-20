\name{MAP}
\alias{MAP}
\title{Velicer's minimum average partial (MAP) test for the number of factors}
\description{Velicer's minimum average partial (MAP) test for determining the number 
	of factors focuses on the common variance in a correlation matrix. It involves 
	a complete principal components analysis followed by the examination of a series 
	of matrices of partial correlations. Specifically, on the first step, the first 
	principal component is partialled out of the correlations between the variables 
	of interest, and the average squared coefficient in the off-diagonals of the 
	resulting partial correlation matrix is computed. On the second step, the first 
	two principal components are partialled out of the original correlation matrix 
	and the average squared partial correlation is again computed. These computations 
	are conducted for k (the number of variables) minus one steps. The average squared 
	partial correlations from these steps are then lined up, and the number of components 
	is determined by the step number in the analyses that resulted in the lowest average 
	squared partial correlation. The average squared coefficient in the original 
	correlation matrix is also computed, and if this coefficient happens to be lower 
	than the lowest average squared partial correlation, then no components should be 
	extracted from the correlation matrix. Statistically, components are retained as 
	long as the variance in the correlation matrix represents systematic variance. 
	Components are no longer retained when there is proportionately more unsystematic 
	variance than systematic variance (see O'Connor, 2000, p. 397).
	
	The MAP test is often more appropriate for factor analyses than it is for principal 
	components analyses. In Velicer's words, "Component analysis has a variety of purposes. 
	It can be used to find a parsimonious description of the total variance of the variables 
	involved; in this case, the [MAP test] is not applicable. Principal component analysis 
	is frequently used to express the variance shared among variables in a set; that is, 
	it is used as kind of a factor analysis" (1976, p. 321). "... if component analysis is 
	employed as an alternative to factor analysis or as a first-stage solution for factor 
	analysis, the stopping rule proposed here would seem the most appropriate." (1976, p. 326)'	
	}
\usage{MAP(data, corkind, verbose)}
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
 \item{eigenvalues}{eigenvalues}
 \item{avgsqrs}{Velicer's average squared correlations}
 \item{nfMAP}{number of factors according to the original (1976) MAP test}
 \item{nfMAP4}{number of factors according to the revised (2000) MAP test}
}
\references{
 {Velicer, W. F. (1976). Determining the number of components from the matrix of partial 
 	correlations. \emph{Psychometrika, 41,} 321-327.}
\cr\cr {Velicer, W. F., Eaton, C. A., and Fava, J. L. (2000). Construct explication 
	through factor or component analysis: A review and 
 evaluation of alternative procedures for determining the number of factors or 
 components. In R. D. Goffin & E. Helmes, eds., \emph{Problems and solutions in human 
 assessment (p.p. 41-71).} Boston: Kluwer.}
\cr\cr {O'Connor, B. P. (2000). SPSS and SAS programs for determining 
 the number of components using parallel analysis and Velicer's 
 MAP test. \emph{Behavior Research Methods, Instrumentation, and
 Computers, 32,} 396-402.}
}
\author{Brian P. O'Connor }
\examples{
# MAP test on the Harman correlation matrix
MAP(data_Harman, corkind='pearson', verbose=TRUE)

\donttest{
# MAP test on the Rosenberg Self-Esteem Scale (RSE) using Pearson correlations
MAP(data_RSE, corkind='pearson', verbose=TRUE)

# MAP test on the Rosenberg Self-Esteem Scale (RSE) using polychoric correlations
MAP(data_RSE, corkind='polychoric', verbose=TRUE)

# MAP test on the NEO-PI-R data
MAP(data_NEOPIR, verbose=TRUE)
}
}
