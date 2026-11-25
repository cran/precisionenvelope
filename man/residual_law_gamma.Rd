\name{residual_law_gamma}
\alias{residual_law_gamma}
\title{A Heteroskedastic Law for the Residual Variance}
\description{
Returns equally weighted quantiles of a gamma distribution with a given
mean and coefficient of variation, for use as the argument \code{b} of
\code{\link{precision_envelope}} when the residual variance
\eqn{b_\theta(W)} is expected to vary across covariate profiles.
}
\usage{
residual_law_gamma(cv, n = 1000L, mean = 1)
}
\arguments{
  \item{cv}{Coefficient of variation of \eqn{b_\theta(W)}, non-negative.
    Zero gives a constant residual variance.}
  \item{n}{Number of quantiles.}
  \item{mean}{Mean of the returned values.}
}
\details{
The values are the gamma quantiles at \eqn{(i - 0.5)/n}, \eqn{i = 1, \dots, n},
rescaled to have exactly the requested mean. The gamma family is a
convenient planning assumption, not a model the method requires; any vector
of plausible residual variances, for example fitted values from a pilot
study, can be passed to \code{\link{precision_envelope}} instead.
}
\value{
A numeric vector of length \code{n}.
}
\seealso{
\code{\link{precision_envelope}}
}
\examples{
b <- residual_law_gamma(cv = 0.5)
c(mean = mean(b), cv = sd(b) / mean(b))
}
\keyword{distribution}
