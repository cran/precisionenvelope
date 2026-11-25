\name{envelope_bootstrap}
\alias{envelope_bootstrap}
\alias{print.envelope_band}
\alias{as.data.frame.envelope_band}
\title{Simultaneous Refitted-Bootstrap Band and Certified Phases}
\description{
Resamples units, refits every nuisance function and recomputes the
guaranteed and possible margins at each grid rate, and forms a
simultaneous band from the maximum studentized deviation. A rate is
certified as guaranteed, allocation-sensitive or inadmissible only when the
band excludes the other phases; otherwise it is reported as indeterminate.
The band is also inverted into confidence intervals for the threshold
missingness rates.
}
\usage{
envelope_bootstrap(fit, R = 999L, alpha = 0.05,
                   failure = c("abort", "discard"),
                   min_success_fraction = 0.75, min_success_count = 40L)

\method{print}{envelope_band}(x, digits = 3, \dots)

\method{as.data.frame}{envelope_band}(x, row.names = NULL, optional = FALSE, \dots)
}
\arguments{
  \item{fit}{An object returned by \code{\link{envelope_fit}} with a
    requirement.}
  \item{R}{Number of bootstrap refits.}
  \item{alpha}{One minus the simultaneous confidence level.}
  \item{failure}{\code{"abort"} (default) stops unless every refit
    succeeds. \code{"discard"} drops failed refits, provided that at least
    \code{min_success_count} refits and a fraction \code{min_success_fraction}
    of them succeed; the band is then conditional on success.}
  \item{min_success_fraction, min_success_count}{Thresholds used with
    \code{failure = "discard"}.}
  \item{x}{An object returned by \code{envelope_bootstrap}.}
  \item{row.names, optional}{Not used; present for compatibility with the
    generic.}
  \item{digits}{Number of decimal places to print.}
  \item{\dots}{Ignored.}
}
\details{
With \eqn{S} successful refits, the critical value is the order statistic of
rank \eqn{\lceil (S + 1)(1 - \alpha) \rceil} of the maximum absolute
studentized deviation, a conservative finite-sample choice; at least
\eqn{\lceil (1 - \alpha)/\alpha \rceil} successful refits are needed. The
guaranteed and possible margins are nonincreasing in the missingness rate,
and the band is tightened under this constraint before it is inverted into
intervals for \eqn{q_G^*} and \eqn{q_P^*} that respect
\eqn{q_G^* \le q_P^*}.

Each refit repeats the full cross-fitted estimation, so the computation time
is about \code{R} times that of \code{\link{envelope_fit}}. Call
\code{\link{set.seed}} beforehand for reproducible output.
}
\value{
An object of class \code{"envelope_band"}, a list with components
\item{band}{A data frame with one row per grid rate: \code{p}, \code{q},
  the estimated margins with their simultaneous lower and upper limits, and
  the certified \code{phase}.}
\item{critical_value, rank, alpha}{The critical value, its rank and the
  level.}
\item{attempted, successful, failure}{Counts of refits and the failure
  policy.}
\item{feasibility}{For an absolute requirement, whether it is certified to
  hold with complete data; otherwise \code{"not applicable"}.}
\item{thresholds}{A matrix with rows \code{guaranteed} and \code{possible}
  and columns \code{lower} and \code{upper}, or \code{NULL} when phases
  cannot be certified.}
\item{fit}{The input fit.}
}
\references{
Nakamura, S. (2026). Sharp estimand-specific precision envelopes for
blockwise missing data under missing at random. \emph{Stats}, \bold{9}(5),
110. \doi{10.3390/stats9050110}
}
\seealso{
\code{\link{envelope_fit}}
}
\examples{
set.seed(2)
n <- 400
w1 <- rnorm(n)
y <- 10 + 2 * w1 + rnorm(n, sd = exp(0.3 * w1))
seen <- rbinom(n, 1, plogis(1 + 0.8 * w1))
dat <- data.frame(y = ifelse(seen == 1, y, NA), w1 = w1)
fit <- envelope_fit(dat, outcome = "y", covariates = "w1",
                    p = c(0.6, 0.8, 1.0), lower = 0.1, target = 0.8)

# A small number of refits keeps the example fast; use R = 999 in practice.
band <- envelope_bootstrap(fit, R = 39)
band
\donttest{
band_full <- envelope_bootstrap(fit, R = 999)
as.data.frame(band_full)
}
}
\keyword{models}
