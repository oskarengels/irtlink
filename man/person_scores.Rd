\name{person_scores}
\alias{person_scores}
\title{Person parameter estimates for an em calibration}
\usage{
person_scores(calib, method = c("eap", "wle"))
}
\arguments{
\item{calib}{An \code{irtlink_calib} object fitted with \code{engine = "em"}.}

\item{method}{\code{"eap"} (default) or \code{"wle"}.}
}
\value{
data.frame with \code{group}, \code{person} (row number within the
group's data), \code{est}, and \code{se}.
}
\description{
Computes person scores from a fitted \code{\link{calibrate}} object of
the internal EM engine (\code{engine = "em"}): \code{"eap"} gives the
expected a posteriori (EAP) score with its posterior standard
deviation under the group's
fitted latent distribution, \code{"wle"} the weighted likelihood
estimate (WLE, Warm, 1989) with the standard error
\eqn{1/\sqrt{I(\hat{\theta})}}, where \eqn{I(\theta)} is the test
information over the person's observed items. Supported for
\code{calibration = "separate"} and \code{"concurrent"} with stored
models (\code{keep_models = TRUE}). Polytomous calibrations (GPCM,
PCM) support \code{"eap"}.
}
\examples{
set.seed(1)
d <- sim_trend_data(n_groups = 2, N = 200, I = 8, overlap = 1,
                    dif = "none")
cal <- calibrate(d)
head(person_scores(cal, method = "eap"))
head(person_scores(cal, method = "wle"))
}
