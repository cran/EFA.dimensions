\name{INTERNAL.CONSISTENCY}
\alias{INTERNAL.CONSISTENCY}
\title{Internal consistency reliability coefficients}
\description{Internal consistency reliability coefficients}
\usage{INTERNAL.CONSISTENCY(data, extraction = 'ML', reverse_these = NULL, 
	                        auto_reverse = TRUE, verbose=TRUE, factormodel)}  
\arguments{
  \item{data}{
  \code{}An all-numeric dataframe where the rows are cases & the columns are the variables.}

  \item{extraction}{ 
  \code{}(optional) The factor extraction method to be used in the omega computations. The options are:
         'ML' for maximum likelihood (the default); and
  	     'PAF' for principal axis / common factor analysis.}
                  
  \item{reverse_these}{
  \code{}(optional) A vector of the names of items that should be reverse-coded}

  \item{auto_reverse}{
  \code{}(optional) Should reverse-coding of items be conducted when warranted? TRUE (default) or FALSE}

  \item{verbose}{
  \code{}(optional) Should detailed results be displayed in console? TRUE (default) or FALSE}
  
  \item{factormodel}{ 
  \code{}(Deprecated.) Use 'extraction' instead.}
}
\details{
When 'auto_reverse = TRUE', the item loadings on the first principal component
are computed and items with negative loadings are reverse-coded.

If error messages are produced, try using 'auto_reverse = FALSE'.
	
If item names are provided for the 'reverse_these' argument, then auto_reverse is not conducted.
	
The following helpful descriptions of Cronbach's alpha and of omega total 
are direct quotes from McNeish (2017, pp. 414-417):

\strong{Cronbach's Alpha}

"One can interpret the value of Cronbach's alpha in one of many
different ways:

1. Cronbach's alpha is the correlation of the scale of interest with
another scale of the same length that intends to measure the same
construct, with different items, taken from the same hypothetical pool
of items (Kline, 1986).

2. The square root of Cronbach's alpha is an estimate of the
correlation between observed scores and true scores (Nunnally &
Bernstein, 1994).

3. Cronbach's alpha is the proportion of the variance of the scale
that can be attributed to a common source (DeVellis, 1991).

4. Cronbach's alpha is the average of all possible split-half
reliabilities from the set of items (Pedhazur & Schmelkin, 1991). (It
is important to note the correlation between the two parts is not the
split half reliability, but is used to find the split half reliability
found by the Spearman-Brown prophecy formula.)

Under certain assumptions, Cronbach's alpha is a consistent estimate
of the population internal consistency; however, these assumptions are
quite rigid and are precisely why methodologists have argued against
the use of Cronbach's alpha.

The assumptions of Cronbach's alpha are:

1. The scale adheres to tau equivalence, i.e., that each item on a
scale contributes equally to the total scale score. Tau equivalence
tends to be unlikely for most scales that are used in empirical
research  some items strongly relate to the construct while some are
more weakly related.

2. Scale items are on a continuous scale and normally distributed.
Cronbach's alpha is largely based on the observed covariances (or
correlations) between items. In most software implementations of
Cronbach's alpha (such as in SAS and SPSS), these item covariances are
calculated using a Pearson covariance matrix. A well-known assumption
of Pearson covariance matrices is that all variables are continuous in
nature. Otherwise, the elements of the matrix can be substantially
biased downward. However, it is particularly common for psychological
scales to contain items that are discrete (e.g., Likert or binary
response scales), which violates this assumption. If discrete items
are treated as continuous, the covariance estimates will be
attenuated, which ultimately results in underestimation of Cronbach's
alpha because the relations between items will appear smaller than
they actually are. To accommodate items that are not on a continuous
scale, the covariances between items can instead be estimated with a
polychoric covariance (or correlation) matrix rather than with a
Pearson covariance matrix. Polychoric covariance matrices assume that
there is an underlying normal distribution to discrete responses.

3. The errors of the items do not covary. Correlated errors occur when
sources other than the construct being measured cause item responses
to be related to one another.

4. The scale is unidimensional. Though Cronbach's alpha is sometimes
thought to be a measure of unidimensionality because its colloquial
definition is that it measures how well items stick together,
unidimensionality is an assumption that needs to be verified prior to
calculating Cronbach's alpha rather than being the focus of what
Cronbach's alpha measures. Internal consistency is necessary for
unidimensionality but that internal consistency is not sufficient for
demonstrating unidimensionality. That is, items that measure different
things can still have a high degree of interrelatedness, so a large
Cronbach's alpha value does not necessarily guarantee that the scale
measures a single construct. As a result, violations of
unidimensionality do not necessarily bias estimates of Cronbach's
alpha. In the presence of a multidimensional scale, Cronbach's alpha
may still estimate the interrelatedness of the items accurately and
the interrelatedness of multidimensional items can in fact be quite
high."

\strong{Omega total}

"Omega total is an internal consistency coefficient that assumes that
the scale is unidimensional. Omega estimates the reliability for the
composite of items on the scale (which is conceptually similar to
Cronbach's alpha). Under the assumption that the construct variance is
constrained to 1 and that there are no error covariances, omega total
is calculated from factor analysis output (loadings and
error/uniqueness values). Tau equivalence is no longer assumed and the
potentially differential contribution of each item to the scale must
be assessed. Omega total is a more general version of Cronbach's alpha
and actually subsumes Cronbach's alpha as a special case. More simply,
if tau equivalence is met, omega total will yield the same result as
Cronbach's alpha but omega total has the flexibility to accommodate
congeneric scales, unlike Cronbach's alpha."

\strong{Root Mean Square Residual (rmsr)}

rmsr is an index of the overall badness-of-fit.
It is the square root of the mean of the squared residuals (the residuals being
the simple differences between original correlations and the correlations
implied by the N-factor model). rmsr is 0 when there is perfect model fit. A
value less than .08 is generally considered a good fit. The rmsr coefficient is 
included in the internal consistency output as an index of the degree of fit of 
a one-factor model to the item data. 

\strong{Standardized Cronbach's Alpha}

Standardized alpha should be used when items have different scale
ranges, e.g., some items are 1-to-7, and other items are 1-to-4, or
1-to-100. Regular alpha is based on covariances, whereas standardized
alpha is based on correlations, wherein the items have identical
standard deviations. Items in different metrics should be standardized
before computing scale scores.
}
\value{A list with the following elements:
 \item{int.consist_scale}{A vector with the scale omega, Cronbach's alpha, 
 	standardized Cronbach's alpha, the mean of the off-diagonal correlations, 
 	the median of the off-diagonal correlations, and
 	the rmsr fit coefficient for a 1-factor model}
 \item{int.consist_dropped}{A matrix of the int.consist_scale values for when each item, 
 	in turn, is int.consist_dropped from the analyses}
 \item{item_stats}{The item means, standard deviations, and item-total correlations}
 \item{resp_opt_freqs}{The response option frequencies}
 \item{resp_opt_props}{The response option proportions}
 \item{new_data}{The data that was used for the analyses, including any item reverse-codings}
}
\references{ 
       {Flora, D. B. (2020). Your coefficient alpha is probably wrong, but which 
       	coefficient omega is right? A tutorial on using R to obtain better reliability estimates. 
       	\emph{Advances in Methods and Practices in Psychological Science, 3(4),} 484501.}
\cr\cr {McNeish, D. (2018). Thanks coefficient alpha, we'll take it from here. 
	    \emph{Psychological Methods, 23(3),} 412433.}
\cr\cr {Revelle, W., & Condon, D. M. (2019). Reliability from alpha to omega: A tutorial. 
	    \emph{Psychological Assessment, 31(12),} 13951411.}
}
\author{Brian P. O'Connor }
\examples{
# Rosenberg Self-Esteem scale items -- without reverse-coding
INTERNAL.CONSISTENCY(data_RSE, extraction = 'PAF', 
                     reverse_these = NULL, auto_reverse = FALSE, verbose=TRUE)

# Rosenberg Self-Esteem scale items -- with auto_reverse-coding
INTERNAL.CONSISTENCY(data_RSE, extraction = 'PAF',
                     reverse_these = NULL, auto_reverse = TRUE, verbose=TRUE)

# Rosenberg Self-Esteem scale items -- another way of reverse-coding
INTERNAL.CONSISTENCY(data_RSE, extraction = 'PAF',
                     reverse_these = c('Q1','Q2','Q4','Q6','Q7'), verbose=TRUE)
}
