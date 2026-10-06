\name{dependent_vcov}
\alias{dependent_vcov}
\title{Item parameter covariance for dependent person samples}
\usage{
dependent_vcov(x, ids, method = c("stacked", "bootstrap"), ...)
}
\arguments{
\item{x}{For \code{method = "stacked"} an \code{irtlink_calib} object
from \code{\link[=calibrate]{calibrate()}} with \code{calibration = "separate"} and
\code{keep_models = TRUE}. For \code{method = "bootstrap"} the list of
response matrices as passed to \code{\link[=calibrate]{calibrate()}}.}

\item{ids}{List with one vector of person identifiers per group, in
the row order of the response data. Identifiers shared between
groups mark the same person. Each identifier may occur at most once
per group. \code{\link[=sim_trend_data]{sim_trend_data()}} stores suitable
identifiers in its \code{person_ids} attribute.}

\item{method}{\code{"stacked"} (default) computes the closed-form
stacked sandwich covariance. \code{"bootstrap"} estimates the same
matrix by a nonparametric bootstrap over persons (e.g., Efron,
1979).}

\item{...}{Passed to the chosen estimator. The stacked sandwich
accepts \code{h} (step size of the observed information, default
\code{1e-5}). The bootstrap accepts \code{B} (number of replicates,
default 200), \code{seed}, \code{model}, \code{cores}, and
\code{draws} (previously computed draws).}
}
\value{
A list whose element \code{vcov} holds the joint covariance of the
stacked item parameters, with row names that join group number, item
name, and parameter as in \verb{"g1:I01:a"} and
\verb{"g1:I01:b"}. It can be passed directly as \code{vcov}
to \code{\link[=linking_error]{linking_error()}} or as \code{le_vcov} to
\code{\link[=estimate_trend]{estimate_trend()}}. Both methods also return
\code{vcov_items} (per-item covariance blocks) and \code{overlap}, a
matrix whose entry for two groups counts the persons who appear in
both. The stacked method additionally returns reference matrices
without the dependence (\code{vcov_nocross},
\code{vcov_independence}) and a \code{score_check}. The bootstrap
additionally returns its \code{draws} and records replicates that
were excluded: a replicate is excluded when a calibration did not
converge or when a resampled group answered an item entirely correct
or entirely wrong, which leaves no response variance to calibrate.
}
\description{
Joint covariance matrix of the item parameter estimates of all
separately calibrated groups when the same persons were tested in
more than one group, as in longitudinal designs that follow a
cohort. Persons shared between groups make the item parameter
estimates of these groups dependent, and the cross-group blocks of
this covariance carry exactly that dependence into
\code{\link[=linking_error]{linking_error()}}. The stacked method treats the separate
calibrations as one stacked M-estimation problem with the
block-diagonal observed information as bread and the crossproduct of
the person-level scores as meat. The bootstrap method redraws
persons with their complete longitudinal records and recalibrates
every group per replicate.
}
\examples{
\donttest{
d <- sim_trend_data(n_groups = 2, N = 300, I = 10,
                    person_overlap = 0.5)
cal <- calibrate(d, keep_models = TRUE)
dv <- dependent_vcov(cal, attr(d, "person_ids"))
dv$overlap
}
}
\references{
Efron, B. (1979). Bootstrap methods: Another look at the jackknife.
The Annals of Statistics, 7(1), 1-26. \doi{10.1214/aos/1176344552}
}
