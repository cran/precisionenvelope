\name{recovery_plan}
\alias{recovery_plan}
\alias{print.recovery_plan}
\title{Minimum-Cost Targeted Recovery of Missing Outcomes}
\description{
Finds the least costly increase of the observation probabilities, from a
current observation law, that brings the residual component
\eqn{E\{b_\theta(W)/\pi(W)\}} down to a required value. The solution is the
cost-adjusted square-root rule
\eqn{\pi(W) = \min[1, \max\{e_0(W), \{\eta\, b_\theta(W)/c(W)\}^{1/2}\}]},
where \eqn{e_0} is the current observation probability and \eqn{c} the cost
of raising it by one unit.
}
\usage{
recovery_plan(b, baseline, A = NULL, target = NULL, variance_limit = NULL,
              residual_target = NULL, cost = 1, weights = NULL)

\method{print}{recovery_plan}(x, digits = 3, \dots)
}
\arguments{
  \item{b}{Residual variances \eqn{b_\theta(W)} of the units or strata.}
  \item{baseline}{Current observation probabilities \eqn{e_0(W)} in
    \eqn{(0, 1]}, one per entry of \code{b}.}
  \item{A}{Between-covariate variance component, needed with \code{target}
    or \code{variance_limit}.}
  \item{target}{Retention requirement in \eqn{(0, 1]}.}
  \item{variance_limit}{Absolute limit on the observed-data efficient
    variance.}
  \item{residual_target}{Required value of the residual component itself.
    Supply exactly one of \code{target}, \code{variance_limit} and
    \code{residual_target}.}
  \item{cost}{Positive marginal cost of raising the observation probability
    by one unit, a single value or one value per entry of \code{b}.}
  \item{weights}{Optional masses of the entries of \code{b}.}
  \item{x}{An object returned by \code{recovery_plan}.}
  \item{digits}{Number of decimal places to print.}
  \item{\dots}{Ignored.}
}
\details{
For a unit that is currently missing, the planned recovery probability is
\eqn{\{\pi(W) - e_0(W)\}/\{1 - e_0(W)\}}. The plan treats \eqn{b_\theta},
\eqn{e_0} and the costs as known; with estimated inputs it is a plug-in
planning calculation and carries no inferential guarantee for an adaptive
recovery design.
}
\value{
An object of class \code{"recovery_plan"}, a list with components
\item{status}{\code{"planned"}, \code{"already meets the requirement"},
  \code{"requires every informative outcome"} or \code{"infeasible"} (the
  requirement cannot be met even after recovering every outcome).}
\item{units}{A data frame with \code{b}, \code{baseline}, \code{cost}, the
  planned probability \code{final}, the recovery probability among missing
  units \code{recovery_among_missing} and the initial priority
  \eqn{b/(c\, e_0^2)}.}
\item{baseline_rate, final_rate, added_rate, added_cost}{Mean observation
  probabilities before and after recovery, their difference and the expected
  added cost per unit.}
\item{C_complete, C_baseline, C_target, C_achieved}{The residual component
  with complete data, at baseline, as required and as achieved.}
\item{eta}{The multiplier of the square-root rule.}
}
\references{
Nakamura, S. (2026). Sharp estimand-specific precision envelopes for
blockwise missing data under missing at random. \emph{Stats}, \bold{9}(5),
110. \doi{10.3390/stats9050110}
}
\seealso{
\code{\link{precision_envelope}}
}
\examples{
# Two equally sized strata; the high-variance stratum is observed half the
# time. Recover enough of it to reach 80\% retention.
plan <- recovery_plan(b = c(0.1, 1.9), baseline = c(1, 0.5), A = 1,
                      target = 0.80)
plan
plan$units
}
\keyword{optimize}
