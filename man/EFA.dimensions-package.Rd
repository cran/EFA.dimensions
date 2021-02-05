\name{EFA.dimensions-package}

\alias{EFA.dimensions-package}

\docType{package}

\title{EFA.dimensions}

\description{This package provides exploratory factor analysis-related functions for 
	assessing dimensionality.
	
	There are 11 functions for determining the number of 
	factors (DIMTESTS, EMPKC, HULL, MAP, NEVALSGT1, PARALLEL, RAWPAR, ROOTFIT, SALIENT,  
	SCREE_PLOT, SESCREE, and SMT). 
	
	There are functions for conducting principal components analysis (PCA), principal 
	axis factor analysis (PA_FA), maximum likelihood factor analysis (MAXLIKE_FA), 
	image factor analysis (IMAGE_FA), and extension factor analysis (EXTENSION_FA), 
	
	Varimax rotation (VARIMAX), promax rotation (PROMAX), and Procrustes rotations 
	(PROCRUSTES) can be performed. 
	
	The analyses can be conducted using raw data or correlation matrices as input.
	
	The analyses can be conducted using Pearson correlations, Kendall correlations, 
	Spearman correlations, Goodman-Kruskal gamma correlations (Thompson, 2006), 
	or polychoric correlations (using the psych and polychor packages). 
	
	Additional functions focus on the factorability of a correlation matrix (FACTORABILITY),  
	the congruences between factors from different datasets (CONGRUENCE), the assessment 
	of local independence (LOCALDEP), and on the assessment of factor solution 
	complexity (COMPLEXITY). 
}
\references{
{Auerswald, M., & Moshagen, M. (2019). How to determine the number of factors to 
	retain in exploratory factor analysis: A comparison of extraction methods under 
	realistic conditions. \emph{Psychological Methods, 24(4),} 468-491.}
\cr\cr {Fabrigar, L. R., & Wegener, D. T. (2012). \emph{Exploratory factor analysis.} 
  New York, NY: Oxford UNiversity Press. ISBN:978-0-19-973417-7}
\cr\cr {Field, A., Miles, J., & Field, Z. (2012). \emph{Discovering statistics using R.}
  Los Angeles, CA: Sage. ISBN:978-1-4462-0045-2}
\cr\cr {O'Connor, B. P. (2000). SPSS and SAS programs for determining 
  the number of components using parallel analysis and Velicer's 
  MAP test. \emph{Behavior Research Methods, Instrumentation, and
  Computers, 32,} 396-402.}
\cr\cr {O'Connor, B. P. (2001). EXTENSION: SAS, SPSS, and MATLAB programs for extension analysis. 
  \emph{Applied Psychological Measurement, 25,} p. 88. doi:10.1177/01466216010251011.}
\cr\cr {Thompson, L. A. 2007. R (and S-PLUS) Manual to Accompany Agresti's Categorical Data 
	Analysis (2002) 2nd edition.}
}
