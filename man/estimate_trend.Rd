\name{estimate_trend}
\alias{estimate_trend}
\title{Estimate a longitudinal IRT trend}
\usage{
estimate_trend(
  data_list,
  calibration = c("separate", "concurrent", "fixed", "regularized"),
  approach = c("chain", "joint", "joint_restricted"),
  link = c("haberman", "haebara", "sl", "mm", "mgm"),
  variant = c("simultaneous", "pairwise"),
  pow = 2,
  weights = "uniform",
  use_intercepts = TRUE,
  item_weights = c("inverse_admin", "uniform"),
  dif = c("none", "purify"),
  le = c("none", "jackknife", "jackknife_bc", "ajk", "ajk_bc", "sandwich_esw",
    "sandwich_osw", "sandwich_bosw", "sandwich_esw_bc", "sandwich_osw_bc",
    "sandwich_bosw_bc", "sandwich_osw_analytic"),
  le_vcov = NULL,
  le_h = 1e-05,
  ref = NULL,
  model = c("2PL", "1PL", "GPCM", "PCM"),
  eps = 0.001,
  ...
)
}
\arguments{
\item{data_list}{List of per-group data.frames or matrices with
dichotomous responses.}

\item{calibration}{How the groups are calibrated. \code{"separate"}
fits every group on its own metric and links afterwards.
\code{"concurrent"} fits one multiple group model under full or
partial invariance. \code{"fixed"} is fixed item parameter
calibration (FIPC). \code{"regularized"} penalizes group differences
of the item difficulties with a smooth information criterion.}

\item{approach}{How the groups are composed: \code{"chain"} over
successive pairs, \code{"joint"} in one fit over all groups, or
\code{"joint_restricted"} over the growing window of the first t
groups. The approach applies to \code{"separate"} and
\code{"concurrent"} calibration alike. It is ignored for
\code{"fixed"}, which is sequential by construction.}

\item{link}{Linking method for \code{calibration = "separate"}:
\code{"haberman"}, \code{"haebara"}, \code{"sl"}, \code{"mm"}, or
\code{"mgm"}. \code{"mm"} and \code{"mgm"} need
\code{approach = "chain"}.}

\item{variant}{\code{"simultaneous"} (default) or \code{"pairwise"}
form of the joint linking methods, as in
\code{\link[=link_groups]{link_groups()}}.}

\item{pow, weights, use_intercepts, item_weights}{Method options
passed to the linking step (see \code{\link[=link_groups]{link_groups()}}).}

\item{dif}{\code{"none"} (default) keeps all common items in the
link. \code{"purify"} first runs the iterative purification of
\code{\link[=detect_dif]{detect_dif()}}: items flagged as showing differential item
functioning (DIF) or item parameter drift (IPD) leave the
link, and the trend is estimated on the remaining items.}

\item{le}{Linking error estimator passed to
\code{\link[=linking_error]{linking_error()}}, or \code{"none"} (default). The
\verb{*_bc} methods need an item parameter covariance (see
\code{le_vcov}).}

\item{le_vcov}{Item parameter covariance for the \verb{*_bc} methods.
When omitted, the per-group covariances of the calibration are used,
which treat the person samples as independent. For dependent person
samples pass the covariance from
\code{\link[=dependent_vcov]{dependent_vcov()}}.}

\item{le_h}{Step size passed to \code{\link[=linking_error]{linking_error()}}.}

\item{ref}{Reference group whose metric carries the reported trend
(\code{mu = 0}, \code{sigma = 1} for this group). Default \code{NULL}
uses the first group. All linking error estimators report their
errors in the reference metric.}

\item{model}{\code{"2PL"} (default), \code{"1PL"}, \code{"GPCM"}, or
\code{"PCM"}. The polytomous models run the separate path with
\code{dif = "none"}. The 3PL is currently available in
\code{\link[=calibrate]{calibrate()}} with separate calibration and
not yet in this convenience wrapper.}

\item{eps}{Tuning constant of the regularized calibration.}

\item{...}{Passed to \code{\link[=calibrate]{calibrate()}} and
\code{\link[=detect_dif]{detect_dif()}}, e.g. person sampling
weights as \code{pweights}.}
}
\value{
An object of class \code{irtlink_trend}: a list with \code{trend}
(data.frame \verb{group, mu, sigma, se, le, te} and, for the
bias-corrected methods, \verb{le_bc, te_bc}), \code{calibration},
\code{approach}, \code{link}, \code{calib}, \code{dif},
\code{link_obj}, \code{le_obj}, and \code{call}.
}
\description{
Runs the complete pipeline in one call: calibrates the groups,
optionally detects and removes DIF/IPD items, links the groups, and
attaches the requested error estimates. \code{calibration} chooses
the calibration, \code{approach} the composition of the groups
(also for concurrent calibration), and \code{link} chooses the
transformation when the groups were calibrated separately. With a
bias-corrected linking error method the trend table decomposes the
uncertainty into standard error (se), linking error (le and le_bc),
and total error (te and te_bc).
}
\examples{
set.seed(1)
d <- sim_trend_data(n_groups = 2, N = 200, I = 10, overlap = 1,
                    dif = "none", mu = c(0, 0.3), sigma = c(1, 1))
tr <- estimate_trend(d, calibration = "separate", approach = "chain",
                     link = "mgm", le = "jackknife")
tr$trend
}
