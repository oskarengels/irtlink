\name{calibrate}
\alias{calibrate}
\title{Calibrate IRT models per group}
\usage{
calibrate(
  data_list,
  model = c("2PL", "1PL", "3PL", "GPCM", "PCM"),
  calibration = c("separate", "concurrent", "fixed", "regularized"),
  engine = "em",
  free_items = NULL,
  eps = 0.001,
  pweights = NULL,
  keep_models = TRUE,
  keep_vcov = FALSE,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{data_list}{List of data.frames or matrices (one per group) with
dichotomous responses (0/1/NA), or integer categories 0, 1, ... for
the polytomous models, and item names as column names.}

\item{model}{\code{"2PL"} (default), \code{"1PL"}, \code{"3PL"},
\code{"GPCM"}, or \code{"PCM"}.
The 3PL estimates the guessing parameter with a beta prior and
requires \code{calibration = "separate"}. The guessing estimates are
carried in the \code{c} column of the canonical \verb{$ipars} table
and enter the Haebara and Stocking-Lord linking criteria.
\code{"GPCM"} is the generalized partial credit model (Muraki, 1992)
for integer response categories 0, 1, ..., and \code{"PCM"} its
partial credit restriction with all discriminations fixed at 1
(Masters, 1982). Both require \code{calibration = "separate"}. Their
\verb{$ipars} rows carry the discrimination \code{a}, the item
location \code{b} (the mean of the step difficulties), and centered
threshold columns \code{tau1, tau2, ...} with \code{b + tau_v} the
step difficulties.}

\item{calibration}{Calibration strategy. \code{"separate"} (default) fits one
model per group and the trend is recovered afterwards by linking. The other
three yield the trend directly: \code{"concurrent"} fits a multiple group model
with invariant item parameters (full, or partial via \code{free_items}),
\code{"fixed"} is sequential fixed-item-parameter calibration (FIPC), and
\code{"regularized"} is SBIC with a smooth-L0 penalty on
differential item functioning (DIF).}

\item{engine}{Estimation backend. The internal EM engine (\code{"em"}) is
the only backend and covers all four calibration strategies, with
control options via \code{control = list(...)} in \code{...}.}

\item{free_items}{Character vector of item names whose difficulty is allowed
to differ across groups (partial invariance). Requires
\code{calibration = "concurrent"}. Default \code{NULL} (full invariance).}

\item{eps}{Smoothing constant(s) for the smooth-L0 BIC penalty in
\code{calibration = "regularized"}. Smaller values yield a sharper L0
approximation. A single value fits once with that constant. A numeric
vector (e.g. \code{c(0.01, 0.001, 0.0001)}) fits a warm-started grid (from
large to small eps) and selects the eps with minimum BIC. Default
\code{0.001}. Ignored for other calibration strategies.}

\item{pweights}{Optional person sampling weights (design weights,
e.g. senate weights): a list with one positive numeric vector per
group, in the row order of \code{data_list}. The weights enter the
likelihood of every calibration strategy, so persons contribute
proportionally to their weight. With \code{keep_vcov = TRUE} the
stored covariances are the inverted weighted information.
\code{\link[=detect_dif]{detect_dif()}} and
\code{\link[=dependent_vcov]{dependent_vcov()}} do not use the
weights yet.}

\item{keep_models}{Keep the fitted model objects? Default \code{TRUE}.
Set \code{FALSE} when memory is a concern (the stored fits carry the
collapsed response patterns).}

\item{keep_vcov}{Compute and store item parameter covariance matrices
in the package's \verb{a,b} parameterization (\verb{a,b,c} for the
3PL). The guessing parameter of the 3PL is estimated with a beta
prior that acts like a penalty, so its covariance describes this
penalized estimator. Currently supported for
\code{calibration = "separate"}. Default \code{FALSE} because the
Hessian computation is costly.}

\item{verbose}{Print estimation progress from the EM engine?
Default \code{FALSE}.}

\item{...}{The only recognized argument is \code{control},
a list of EM options (e.g. \code{theta_nodes}, \code{theta_range},
\code{maxit}, \code{conv}, or \code{lower_a}, which overrides the
default lower bound 0.1 of the discrimination parameter).}
}
\value{
An object of class \code{irtlink_calib}: list with \code{ipars}
(data.frame \verb{group, item, a, b, c}), \code{item_parameters} (data.frame
\verb{group, item, a, b, d, nu, estimated_parameterization} where available),
\code{model}, \code{calibration}, \code{engine}, \code{item_parameterization}, \code{n_groups}, \code{N},
\code{converged}, \code{models}, \code{vcov_ipars}, \code{design}, \code{trend}, \code{free_items},
\code{eps}, and \code{call}.
\verb{$trend} is \code{NULL} for separate calibration and a data.frame
\verb{group, mu, sigma} (group 1 = reference: mu 0, sigma 1) for concurrent,
fixed, and regularized calibration.
\verb{$models} is a per-group list of fitted model objects for
\code{calibration = "separate"}, a single multiple group fit (in a length-1 list)
for \code{calibration = "concurrent"}, and \code{NULL} for \code{"fixed"} or
\code{"regularized"} (or whenever \code{keep_models = FALSE}).
\verb{$vcov_ipars} is a per-group list of covariance matrices for separate
calibrations when \code{keep_vcov = TRUE}, otherwise \code{NULL}.
Row names append \verb{_a} / \verb{_b} (plus \verb{_c} for the 3PL) to the
item name (e.g., \verb{I01_a}), and
\code{\link[=linking_error]{linking_error()}} extracts the \verb{a,b} blocks it needs by name.
\verb{$free_items} is the character vector of items freed under partial
invariance (\code{calibration = "concurrent"}), or
\code{NULL} when no items were freed.
\verb{$eps} is the selected SBIC smoothing constant for
\code{calibration = "regularized"}, or \code{NULL} otherwise.
}
\description{
Fits a unidimensional IRT model for each group on the package's internal
marginal maximum likelihood EM engine (Bock & Aitkin, 1981, with
update formulas following Hanson, 1998) and
collects the item parameters in a unified long format.
The response model is
\deqn{P(X_{pi} = 1 \mid \theta_p) = c_i + (1 - c_i)\,
\mathrm{logit}^{-1}\!\big(a_i(\theta_p - b_i)\big),}
with guessing \eqn{c_i = 0} for the 1PL and 2PL and discrimination
\eqn{a_i = 1} for the 1PL. Separate
calibration fits one model per group. Concurrent calibration fits a
single multiple group model (Hanson & B\enc{é}{e}guin, 2002) and
extracts the latent trend directly from the group distribution
parameters. Fixed item parameter calibration (FIPC, e.g., Kim, 2006) fixes
the group-1 item parameters in the later groups. Regularized
calibration frees a DIF effect per item and group and shrinks it
under a smooth Bayesian information criterion (SBIC, O'Neill &
Burke, 2023, and Robitzsch, 2024).
}
\details{
The strategies constrain the common items differently.
\code{"concurrent"} and \code{"fixed"} hold all common item
parameters equal across groups, so DIF or item parameter drift
(IPD) in a common item biases
their trend estimates. The item bank of FIPC is the first group by
construction, and the reporting metric remains selectable through
the \code{ref} argument of \code{\link[=estimate_trend]{estimate_trend()}}.
\code{"regularized"} frees a difficulty effect per item and group
and shrinks it toward zero, which keeps functionally invariant items
constrained while clearly drifting items escape the penalty. The
discriminations are not regularized. \code{"separate"} places no
constraints during calibration, and DIF/IPD is then handled in the
linking step or by detection and purification
(\code{\link[=detect_dif]{detect_dif()}}).
}
\examples{
\donttest{
d <- sim_trend_data(n_groups = 2, N = 300, I = 10)
cal <- calibrate(d)
head(coef(cal))
}

}
\references{
Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood
estimation of item parameters: Application of an EM algorithm.
Psychometrika, 46(4), 443-459. \doi{10.1007/BF02293801}

Hanson, B. A. (1998). IRT parameter estimation using the EM
algorithm. Technical note.

Hanson, B. A., & B\enc{é}{e}guin, A. A. (2002). Obtaining a common scale for
item response theory item parameters using separate versus concurrent
estimation in the common-item equating design. Applied Psychological
Measurement, 26(1), 3-24. \doi{10.1177/0146621602026001001}

Kim, S. (2006). A comparative study of IRT fixed parameter
calibration methods. Journal of Educational Measurement, 43(4),
355-381. \doi{10.1111/j.1745-3984.2006.00021.x}

Masters, G. N. (1982). A Rasch model for partial credit scoring.
Psychometrika, 47(2), 149-174. \doi{10.1007/BF02296272}

Muraki, E. (1992). A generalized partial credit model: Application
of an EM algorithm. Applied Psychological Measurement, 16(2),
159-176. \doi{10.1177/014662169201600206}

O'Neill, M., & Burke, K. (2023). Variable selection using a smooth
information criterion for distributional regression models.
Statistics and Computing, 33, 71. \doi{10.1007/s11222-023-10204-8}

Robitzsch, A. (2024). Smooth information criterion for regularized
estimation of item response models. Algorithms, 17(4), 153.
\doi{10.3390/a17040153}
}
