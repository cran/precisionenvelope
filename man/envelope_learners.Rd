\name{envelope_learners}
\alias{envelope_learners}
\alias{learner_linear}
\alias{learner_logistic}
\title{Nuisance Learners for Envelope Estimation}
\description{
\code{envelope_learners} collects the three learners used by
\code{\link{envelope_fit}}: one for the outcome mean, one for the residual
variance and one for the observation probability. \code{learner_linear}
and \code{learner_logistic} construct the default least-squares and
logistic regression learners.
}
\usage{
envelope_learners(mean = learner_linear(), variance = learner_linear(),
                  propensity = learner_logistic())

learner_linear()

learner_logistic()
}
\arguments{
  \item{mean}{Learner for the outcome mean \eqn{E(Y \mid W)}, fitted on
    observed outcomes.}
  \item{variance}{Learner for the residual variance, fitted to squared
    residuals.}
  \item{propensity}{Learner for the observation probability, fitted to the
    0/1 observation indicator. It must return probabilities.}
}
\details{
A learner is a function \code{function(x, y, newx)} that fits a model to
the numeric covariate matrix \code{x} and response \code{y} and returns a
numeric vector of predictions for the rows of the matrix \code{newx}. Any
regression method can be wrapped in this form, for example a spline or a
random forest from another package. The defaults standardize the
covariates and fit, with an intercept, ordinary least squares
(\code{learner_linear}) or logistic regression (\code{learner_logistic}).
Fitted residual variances are bounded below and fitted probabilities are
clipped by \code{\link{envelope_fit}}.
}
\value{
\code{envelope_learners} returns a named list of three functions;
\code{learner_linear} and \code{learner_logistic} return a learner.
}
\seealso{
\code{\link{envelope_fit}}
}
\examples{
# A learner that adds squared terms of every covariate.
quadratic <- function(x, y, newx) {
  learner_linear()(cbind(x, x^2), y, cbind(newx, newx^2))
}
learners <- envelope_learners(mean = quadratic, variance = quadratic)
names(learners)
}
\keyword{models}
