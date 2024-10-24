% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{getChiSqStat}
\alias{getChiSqStat}
\title{Compute Chi-Squared Statistics for Variables}
\usage{
getChiSqStat(datX, response)
}
\arguments{
\item{datX}{A matrix or data frame containing predictor variables. It can
consist of both numerical and categorical variables.}

\item{response}{A factor representing the class labels. It must have at least
two levels for the chi-squared test to be applicable.}
}
\value{
A vector of chi-squared statistics, one for each predictor variable
in \code{datX}. For numerical variables, the chi-squared statistic is computed
after binning the variable.
}
\description{
This function calculates the chi-squared statistic for each column of \code{datX}
against the response variable \code{response}. It supports both numerical and
categorical predictors in \code{datX}. For numerical variables, it automatically
discretizes them into factor levels based on standard deviations and mean,
using different splitting criteria depending on the sample size.
}
\details{
For each variable in \code{datX}, the function first checks if the
variable is numerical. If so, it is discretized into factor levels using
either two or three split points, depending on the sample size and the
number of levels in the \code{response}. Missing values are handled by assigning
them to a new factor level.

The chi-squared statistic is then computed between each predictor and the
\code{response}. If the chi-squared test has more than one degree of freedom,
the Wilson-Hilferty transformation is applied to adjust the statistic to a
1-degree-of-freedom chi-squared distribution.
}
\examples{
datX <- data.frame(var1 = rnorm(100), var2 = factor(sample(letters[1:3], 100, replace = TRUE)))
y <- factor(sample(c("A", "B"), 100, replace = TRUE))
getChiSqStat(datX, y)
}
\references{
Loh, W. Y. (2009). Improving the precision of classification
trees. \emph{The Annals of Applied Statistics}, 1710–1737. JSTOR.
}
