\name{precisionenvelope-package}
\alias{precisionenvelope-package}
\alias{precisionenvelope}
\docType{package}
\title{Sharp Precision Envelopes for Outcomes Missing at Random}
\description{
When an outcome, or an all-or-none block of variables, is missing at random
(MAR) given always observed covariates \eqn{W}, the overall missingness rate
does not determine how much precision survives. Two mechanisms with the same
rate can leave very different efficient variances, because the loss depends
on which units go missing. This package computes the range of efficient
variance, and of information retention, that the MAR mechanisms with a given
overall observation rate \eqn{p} and a lower bound \eqn{\pi_L} on the
observation probability can produce. The range is sharp: both endpoints are
attained or approached by admissible mechanisms.

For a scalar target \eqn{\theta} with complete-data efficient influence
function \eqn{\varphi_\theta(Z)}, write
\eqn{a_\theta(W) = E\{\varphi_\theta(Z) \mid W\}} and
\eqn{b_\theta(W) = \mathrm{Var}\{\varphi_\theta(Z) \mid W\}},
\eqn{A_\theta = \mathrm{Var}\{a_\theta(W)\}} and
\eqn{B_\theta = E\{b_\theta(W)\}}. Under a mechanism \eqn{\pi(W)} the
observed-data efficient variance is \eqn{A_\theta + C_\theta(\pi)} with
\eqn{C_\theta(\pi) = E\{b_\theta(W) / \pi(W)\}}, and the information
retention is \eqn{\rho_\theta = (A_\theta + B_\theta) / \{A_\theta +
C_\theta(\pi)\}}. The package reports the smallest and largest
\eqn{C_\theta(\pi)} over the admissible class.

The package has three layers.
\describe{
  \item{Planning}{\code{\link{precision_envelope}} computes the envelope
  from \eqn{A_\theta} (or the share \eqn{R^2} of the influence-function
  variance explained by \eqn{W}) and a law for \eqn{b_\theta(W)}. No data
  are needed. With a precision requirement it classifies each rate as
  guaranteed, allocation-sensitive or inadmissible and returns the threshold
  missingness rates. \code{\link{residual_law_gamma}} supplies a convenient
  law for a heteroskedastic \eqn{b_\theta(W)}.}
  \item{Recovery}{\code{\link{recovery_plan}} gives the minimum-cost
  targeted recovery of missing outcomes that meets a requirement.}
  \item{Estimation}{\code{\link{envelope_fit}} estimates the envelope from
  data on a partially observed continuous outcome, for the mean of that
  outcome, and \code{\link{envelope_bootstrap}} adds a simultaneous
  refitted-bootstrap band over the grid of rates, certified phases and
  confidence intervals for the threshold rates.}
}

The package does not impute missing values, does not improve the estimate of
\eqn{\theta}, and does not test MAR. It states how much precision can survive
under MAR and a stated lower bound on the observation probability.
}
\references{
Nakamura, S. (2026). Sharp estimand-specific precision envelopes for
blockwise missing data under missing at random. \emph{Stats}, \bold{9}(5),
110. \doi{10.3390/stats9050110}
}
\author{
Satoshi Nakamura \email{s-nakamura@med.mie-u.ac.jp}
}
\keyword{package}
