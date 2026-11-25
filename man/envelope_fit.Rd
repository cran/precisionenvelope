\name{envelope_fit}
\alias{envelope_fit}
\alias{print.envelope_fit}
\alias{plot.envelope_fit}
\alias{as.data.frame.envelope_fit}
\title{Estimate the Precision Envelope from Data}
\description{
Estimates, for the mean of a partially observed continuous outcome, the
variance components, the realized information retention and the sharp
precision envelope over a prespecified grid of observation rates. The
nuisance functions are fitted by cross-fitting and the envelope endpoints by
one-step estimators.
}
\usage{
envelope_fit(data, outcome, covariates, observed = NULL, p, lower,
             target = NULL, variance_limit = NULL, folds = 2L,
             learners = envelope_learners(), variance_floor = NULL,
             propensity_clip = c(0.03, 0.995), min_observed_training = 20L)

\method{print}{envelope_fit}(x, digits = 3, \dots)

\method{plot}{envelope_fit}(x, y, \dots)

\method{as.data.frame}{envelope_fit}(x, row.names = NULL, optional = FALSE, \dots)
}
\arguments{
  \item{data}{A data frame with one row per unit.}
  \item{outcome}{Name of the numeric outcome column. Missing values mark
    missing outcomes unless \code{observed} is given.}
  \item{covariates}{Names of the completely observed numeric covariate
    columns \eqn{W}.}
  \item{observed}{Optional name of a 0/1 (or logical) column marking
    observed outcomes.}
  \item{p}{Strictly increasing grid of overall observation rates at which
    the envelope is estimated, each between \code{lower} and 1.}
  \item{lower}{Lower bound \eqn{\pi_L} on the observation probability,
    strictly between 0 and 1. It is an assumption stated by the analyst and
    cannot be checked from the data.}
  \item{target, variance_limit}{Optional precision requirement, as in
    \code{\link{precision_envelope}}. A requirement is needed for
    \code{\link{envelope_bootstrap}}.}
  \item{folds}{Number of cross-fitting folds, at least 2.}
  \item{learners}{A list of nuisance learners, see
    \code{\link{envelope_learners}}. Elements left out take the defaults.}
  \item{variance_floor}{Positive lower bound applied to the fitted residual
    variance. The default is one thousandth of the variance of the observed
    outcomes.}
  \item{propensity_clip}{Lower and upper bounds applied to the fitted
    observation probabilities.}
  \item{min_observed_training}{Minimum number of observed outcomes required
    in each training fold.}
  \item{x}{An object returned by \code{envelope_fit}.}
  \item{y, row.names, optional}{Not used; present for compatibility with the
    generic.}
  \item{digits}{Number of decimal places to print.}
  \item{\dots}{Further arguments passed to \code{\link[graphics]{plot}}, or
    ignored.}
}
\details{
The target is the outcome mean, whose complete-data influence function is
\eqn{Y - \mu}. Then \eqn{a(W) = m(W) - \mu} and \eqn{b(W) = v(W)}, with
\eqn{m(W) = E(Y \mid W)} and \eqn{v(W) = \mathrm{Var}(Y \mid W)}. Within each
cross-fitting fold the outcome mean is fitted on the observed outcomes of
the training folds, the residual variance is fitted to squared residuals
taken from an inner two-way split of the same observed outcomes, and the
observation probability \eqn{e(W)} is fitted on all training rows. The
one-step estimators of \eqn{\mu}, \eqn{A}, \eqn{B} and the realized
\eqn{C = E\{v(W)/e(W)\}} use the augmented inverse-probability weighted
influence functions; the envelope endpoints at each grid rate use the
fitted favorable allocation and the fitted upper quantile of \eqn{v(W)}
for the adverse endpoint.

The reported Wald interval for the realized retention is marginal. For
simultaneous statements over the grid, and for phases certified at a
stated level, use \code{\link{envelope_bootstrap}}. The fitted share of
observation probabilities below \code{lower} is reported as a diagnostic; it
does not validate the lower bound or MAR.

Results depend on the random split into folds; call
\code{\link{set.seed}} beforehand for reproducible output.
}
\value{
An object of class \code{"envelope_fit"}, a list including
\item{mu}{The estimated outcome mean.}
\item{A, B, C_realized, V_full, V_observed}{Estimated variance components,
  complete-data efficient variance and observed-data efficient variance
  under the realized mechanism.}
\item{rho_realized, rho_se, rho_ci}{Realized information retention, its
  standard error and a marginal 95\% Wald interval.}
\item{table}{A data frame with one row per grid rate: \code{p}, \code{q},
  \code{C_favorable}, \code{C_adverse}, \code{rho_favorable},
  \code{rho_adverse} and, with a requirement, \code{margin_guaranteed},
  \code{margin_possible} and the point classification \code{phase}.}
\item{requirement, theta}{The requirement and the vector of margins used by
  \code{\link{envelope_bootstrap}}.}
\item{n, observed_count, observation_rate, share_below_lower}{Sample
  description and the overlap diagnostic.}
\item{nuisance}{The cross-fitted predictions \code{m}, \code{v} and
  \code{e}.}
The object also keeps the specification, the learners and the data needed
for refitting.
}
\references{
Nakamura, S. (2026). Sharp estimand-specific precision envelopes for
blockwise missing data under missing at random. \emph{Stats}, \bold{9}(5),
110. \doi{10.3390/stats9050110}
}
\seealso{
\code{\link{envelope_bootstrap}}, \code{\link{envelope_learners}},
\code{\link{precision_envelope}}
}
\examples{
set.seed(1)
n <- 500
w1 <- rnorm(n)
w2 <- rbinom(n, 1, 0.4)
y <- 10 + 2 * w1 + 1.5 * w2 + rnorm(n, sd = exp(0.3 * w1))
seen <- rbinom(n, 1, plogis(1 + 0.8 * w1))
dat <- data.frame(y = ifelse(seen == 1, y, NA), w1 = w1, w2 = w2)

fit <- envelope_fit(dat, outcome = "y", covariates = c("w1", "w2"),
                    p = c(0.6, 0.7, 0.8, 0.9), lower = 0.1, target = 0.8)
fit
plot(fit)
}
\keyword{models}
