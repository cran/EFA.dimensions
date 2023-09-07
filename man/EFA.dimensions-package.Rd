\name{EFA.dimensions-package}

\alias{EFA.dimensions-package}

\docType{package}

\title{EFA.dimensions}

\description{This package provides exploratory factor analysis-related functions for 
	assessing dimensionality.
\cr\cr  
    There are 11 functions for determining the number of factors (DIMTESTS, EMPKC, 
	HULL, MAP, NEVALSGT1, PARALLEL, RAWPAR, ROOTFIT, SALIENT, SCREE_PLOT, SESCREE, and SMT).
\cr\cr  
    There is a principal components analysis function (PCA), and an exploratory factor 
    analysis function (EFA) with 9 possible factor extraction methods. 
\cr\cr  
	There are 15 possible factor rotation methods that can be used with PCA and EFA. 
\cr\cr  
	The analyses can be conducted using raw data or correlation matrices as input.
\cr\cr  
	The analyses can be conducted using Pearson correlations, Kendall correlations, 
	Spearman correlations, Goodman-Kruskal gamma correlations (Thompson, 2006), 
	or polychoric correlations (using the psych and polychor packages). 
\cr\cr  
	Additional functions focus on the factorability of a correlation matrix (FACTORABILITY),  
	the congruences between factors from different datasets (CONGRUENCE), the assessment 
	of local independence (LOCALDEP), the assessment of factor solution 
	complexity (COMPLEXITY), and internal consistency (INTERNAL.CONSISTENCY). 
}
\references{
{Auerswald, M., & Moshagen, M. (2019). How to determine the number of factors to 
	retain in exploratory factor analysis: A comparison of extraction methods under 
	realistic conditions. \emph{Psychological Methods, 24(4),} 468-491.}
\cr\cr {Field, A., Miles, J., & Field, Z. (2012). \emph{Discovering statistics using R.}
  Los Angeles, CA: Sage. ISBN:978-1-4462-0045-2}
\cr\cr{Mulaik, S. A. (2010). \emph{Foundations of factor analysis} (2nd ed.). Boca Raton, FL: Chapman 
	and Hall/CRC Press, Taylor & Francis Group.}
\cr\cr {O'Connor, B. P. (2000). SPSS and SAS programs for determining 
  the number of components using parallel analysis and Velicer's 
  MAP test. \emph{Behavior Research Methods, Instrumentation, and
  Computers, 32,} 396-402.}
\cr\cr {O'Connor, B. P. (2000). SPSS and SAS programs for determining 
  the number of components using parallel analysis and Velicer's 
  MAP test. \emph{Behavior Research Methods, Instrumentation, and
  Computers, 32,} 396-402.}
\cr\cr {Sellbom, M., & Tellegen, A. (2019). Factor analysis in psychological assessment research: 
	Common pitfalls and recommendations. 
   \emph{Psychological Assessment, 31(12),} 1428-1441. https://doi.org/10.1037/pas0000623}
\cr\cr {Watts, A. L., Greene, A. L., Ringwald, W., Forbes, M. K., Brandes, C. M., Levin-Aspenson, 
	H. F., & Delawalla, C. (2023). Factor analysis in personality disorders research: Modern issues 
	and illustrations of practical recommendations. 
  \emph{Personality Disorders: Theory, Research, and Treatment, 14(1),} 105-117. 
  https://doi.org/10.1037/per0000581}
}

  
