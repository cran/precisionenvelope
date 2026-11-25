\name{precision_envelope}
\alias{precision_envelope}
\alias{print.precision_envelope}
\alias{plot.precision_envelope}
\alias{as.data.frame.precision_envelope}
\title{Sharp Precision Envelope for a Planned Study}
\description{
Computes, for each overall observation rate \eqn{p}, the smallest and largest
efficient variance that a missing-at-random mechanism with observation
probability bounded below by \code{lower} can produce, together with the
value under missing completely at random (MCAR). With a precision
requirement it classifies each rate and returns the threshold missingness
rates. No data are needed: the inputs describe the target and the
covariates at the planning stage.
}
\usage{
precision_envelope(p, lower, A = NULL, r2 = NULL, b = 1, weights = NULL,
                   covariate = c("continuous", "strata"),
                   target = NULL, variance_limit = NULL)

\method{print}{precision_envelope}(x, digits = 3, \dots)

\method{plot}{precision_envelope}(x, y, n = 201L, \dots)

\method{as.data.frame}{precision_envelope}(x, row.names = NULL, optional = FALSE, \dots)
}
\arguments{
  \item{p}{Numeric vector of overall observation rates, each between
    \code{lower} and 1. The missingness rate is \eqn{q = 1 - p}.}
  \item{lower}{Lower bound \eqn{\pi_L} on the observation probability
    \eqn{\pi(W)}, strictly between 0 and 1. Missingness rates up to
    \eqn{1 - \pi_L} are admissible.}
  \item{A}{The between-covariate variance component
    \eqn{A_\theta = \mathrm{Var}\{a_\theta(W)\}}, on the same scale as
    \code{b}. Supply exactly one of \code{A} and \code{r2}.}
  \item{r2}{The share \eqn{R^2 = A_\theta / (A_\theta + B_\theta)} of the
    complete-data influence-function variance explained by \eqn{W}, in
    \eqn{[0, 1)}. For the mean of an outcome \eqn{Y} this is the usual
    \eqn{R^2} of \eqn{Y} on \eqn{W}.}
  \item{b}{Values of the residual variance \eqn{b_\theta(W)}, describing its
    law. A single value (the default) means that \eqn{b_\theta(W)} is
    constant. A vector is read as a sample of units or a discretized law,
    with masses \code{weights}. Only the shape matters for retention, so
    \code{b} can be on any scale when \code{r2} is used.}
  \item{weights}{Optional non-negative masses for the values in \code{b};
    equal masses by default.}
  \item{covariate}{\code{"continuous"} (default) treats the covariate space
    as having no atoms, so that a mechanism can separate units with the same
    value of \eqn{b_\theta(W)}; the adverse endpoint then has the closed
    upper-tail form. \code{"strata"} treats each entry of \code{b} as an
    indivisible stratum of a discrete \eqn{W}, within which every unit
    shares one observation probability; the adverse endpoint is then found
    by exact enumeration, for at most 15 strata.}
  \item{target}{Optional retention requirement \eqn{\rho_0} in
    \eqn{(0, 1]}: the requirement holds when
    \eqn{\rho_\theta \ge \rho_0}.}
  \item{variance_limit}{Optional absolute requirement: the observed-data
    efficient variance must not exceed this value. Supply at most one of
    \code{target} and \code{variance_limit}.}
  \item{x}{An object returned by \code{precision_envelope}.}
  \item{y}{Not used; present for compatibility with the generic.}
  \item{row.names, optional}{Not used; present for compatibility with the generic.}
  \item{digits}{Number of decimal places to print.}
  \item{n}{Number of missingness rates at which the curves are drawn.}
  \item{\dots}{Further arguments passed to \code{\link[graphics]{plot}}, or
    ignored.}
}
\details{
Write \eqn{B_\theta = E\{b_\theta(W)\}}. Over observation laws with
\eqn{\pi_L \le \pi(W) \le 1} and \eqn{E\{\pi(W)\} = p}, the favorable
endpoint \eqn{C^-(p)} minimizes \eqn{E\{b_\theta(W)/\pi(W)\}} and is
attained by the clipped square-root rule
\eqn{\pi(W) = \min[1, \max\{\pi_L, \{b_\theta(W)/\lambda\}^{1/2}\}]}. On a
covariate space without atoms the adverse endpoint is
\deqn{C^+(p) = B_\theta + (1/\pi_L - 1) \int_{1-a}^{1} Q_b(u)\,du,
\qquad a = (1 - p)/(1 - \pi_L),}{C+(p) = B + (1/piL - 1) * integral of the upper a-tail of b, a = (1 - p)/(1 - piL),}
where \eqn{Q_b} is the quantile function of \eqn{b_\theta(W)}. The MCAR
value is \eqn{B_\theta / p}, which always lies between the endpoints.

Information retention is \eqn{\rho = (A_\theta + B_\theta)/(A_\theta + C)},
and the standard error is inflated by the factor \eqn{\rho^{-1/2}}. The
MCAR-equivalent missingness rate of a value \eqn{C} is
\eqn{1 - B_\theta / C}, the MCAR rate that would give the same retention.

A requirement translates into a residual budget \eqn{C_{req}}, equal to
\eqn{(A_\theta + B_\theta)/\rho_0 - A_\theta} for a retention requirement
and to \code{variance_limit} \eqn{- A_\theta} for an absolute one. A rate is
guaranteed when \eqn{C^+(p) \le C_{req}}, so that every admissible mechanism
meets the requirement; allocation-sensitive when only
\eqn{C^-(p) \le C_{req}}, so that some mechanisms meet it and others do not;
and inadmissible when \eqn{C^-(p) > C_{req}}. The thresholds \eqn{q_G^*},
\eqn{q_{MCAR}^*} and \eqn{q_P^*} are the missingness rates at which the
adverse, MCAR and favorable margins reach zero. An absolute requirement that
fails even with complete data is reported as infeasible. A requirement on
the half-width \eqn{h} of a 95\% interval from \eqn{n} units corresponds to
\code{variance_limit = n * (h / qnorm(0.975))^2}.
}
\value{
An object of class \code{"precision_envelope"}, a list with components
\item{table}{A data frame with one row per rate: \code{p}, \code{q}, the
  residual components \code{C_favorable}, \code{C_mcar} and
  \code{C_adverse}, the retentions \code{rho_favorable}, \code{rho_mcar}
  and \code{rho_adverse}, the relative standard error increases
  \code{se_inflation_*}, the MCAR-equivalent missingness rates
  \code{q_equivalent_favorable} and \code{q_equivalent_adverse}, and, when a
  requirement is given, \code{margin_guaranteed}
  (\eqn{C_{req} - C^+}), \code{margin_possible} (\eqn{C_{req} - C^-}) and
  \code{phase}.}
\item{A, B, V_full}{The variance components and the complete-data
  efficient variance \eqn{A_\theta + B_\theta}.}
\item{lower, q_max}{The lower bound and the largest admissible missingness
  rate.}
\item{requirement}{\code{NULL}, or a list with the requirement type, its
  value, \code{C_required}, the retention it implies and whether it can hold
  with complete data.}
\item{thresholds}{\code{NULL}, or a named vector with the threshold
  missingness rates \code{guaranteed}, \code{mcar} and \code{possible}.
  A value equal to \code{q_max} means the margin stays non-negative over the
  whole admissible range; \code{NA} means the requirement fails with
  complete data.}
\item{covariate, b, weights}{The inputs that describe the law of
  \eqn{b_\theta(W)}.}
The \code{plot} method draws the three retention curves against the
missingness rate and returns their values invisibly.
}
\references{
Nakamura, S. (2026). Sharp estimand-specific precision envelopes for
blockwise missing data under missing at random. \emph{Stats}, \bold{9}(5),
110. \doi{10.3390/stats9050110}
}
\seealso{
\code{\link{residual_law_gamma}}, \code{\link{recovery_plan}},
\code{\link{envelope_fit}}
}
\examples{
# Planning with a constant residual variance: 25\% of outcomes expected to be
# missing, covariates explaining 30\% of the outcome variance, and an
# observation probability of at least 0.10 for every covariate profile.
env <- precision_envelope(p = c(0.95, 0.85, 0.75, 0.65), lower = 0.10,
                          r2 = 0.30, target = 0.80)
env

# The analytic illustration of Nakamura (2026, Section 5.1): residual
# variance 0.1 or 1.9 with equal probability, A = 1 and lower bound 0.2.
ill <- precision_envelope(p = c(0.9, 0.75, 0.6), lower = 0.2, A = 1,
                          b = c(0.1, 1.9), target = 0.80)
ill$thresholds
plot(ill)

# A heteroskedastic law with coefficient of variation 0.8.
het <- precision_envelope(p = 0.75, lower = 0.10, r2 = 0.30,
                          b = residual_law_gamma(cv = 0.8))
as.data.frame(het)[, c("q", "rho_favorable", "rho_mcar", "rho_adverse")]
}
\keyword{models}
