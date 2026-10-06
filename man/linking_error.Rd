\name{linking_error}
\alias{linking_error}
\title{Linking error of a linked trend}
\usage{
linking_error(
  link,
  method = c("jackknife", "jackknife_bc", "ajk", "ajk_bc", "sandwich_esw",
    "sandwich_osw", "sandwich_bosw", "sandwich_esw_bc", "sandwich_osw_bc",
    "sandwich_bosw_bc", "sandwich_osw_analytic"),
  vcov = NULL,
  h = 1e-05,
  jk_scale = c("finite_item", "classical"),
  cluster = NULL,
  ...
)
}
\arguments{
\item{link}{An \code{irtlink_link} object returned by
\code{\link[=link_groups]{link_groups()}}.}

\item{method}{\code{"jackknife"} refits the linking without each link
item in turn. \code{"ajk"} replaces the refits by approximate
one-step updates of the estimating equations (Robitzsch, 2024,
2025a, 2025b). \code{"sandwich_esw"}, \code{"sandwich_osw"}, and
\code{"sandwich_bosw"} are sandwich estimators with the expected,
observed, and finite-item corrected observed meat, and
\code{"sandwich_osw_analytic"} evaluates the observed sandwich with
analytic derivatives. The \code{"_bc"} variants apply the
finite-sample bias correction of Robitzsch (2024, 2025a, 2025b) and
require \code{vcov}.}

\item{vcov}{Item parameter covariance matrix or per-group list of
covariance matrices, required for the \code{"_bc"} methods. Use the
\code{vcov_ipars} of \code{\link[=calibrate]{calibrate()}} with
\code{keep_vcov = TRUE} for independent person samples and
\code{\link[=dependent_vcov]{dependent_vcov()}} for dependent ones.}

\item{h}{Step size of the numeric fallback derivatives of the
bias-corrected jackknife. The default analytic derivatives do not
use it.}

\item{jk_scale}{Scaling of the squared delete-one differences for
\code{method = "jackknife"}. \code{"finite_item"} (default) uses the
factor I/(I-1) and \code{"classical"} the factor (I-1)/I of the
classical jackknife.}

\item{cluster}{Optional named character vector mapping item names to
testlet labels, e.g. \code{c(I01 = "T1", I02 = "T1")}. The deleted
jackknife units are then whole testlets instead of single items,
which keeps the linking error valid when items within a testlet
share stimulus material and drift together (cf. the item-unit
jackknife in PISA, Monseur & Berezner, 2007). Items without an
entry stay their own unit. Available for \code{"jackknife"} (every
linking method) and \code{"ajk"} (chain mean-mean and
mean-geometric mean linking and pairwise Haberman linking).
\code{n_reps} then counts units.}

\item{...}{Options overriding the stored linking controls when the
leave-one-item links are refit for \code{"jackknife"}. Ignored with
a warning by the other methods.}
}
\value{
The input \code{irtlink_link} object with its \code{le} slot filled: a
data.frame with \code{group}, \code{le_mu}, \code{le_sigma},
\code{method}, and \code{n_reps}. The \code{"_bc"} methods add
\verb{se_*}, \verb{lebc_*}, \verb{te_*}, and \verb{tebc_*} columns.
}
\description{
Estimates the linking error, the uncertainty of the linked means and
standard deviations that stems from the sampling of the link items
(e.g., Monseur & Berezner, 2007). The delete-one item jackknife
aggregates the refits without each link item around the full
estimate,
\deqn{LE^2(\hat{\mu}_t) = \frac{I}{I-1} \sum_{i=1}^{I}
\left(\hat{\mu}_{t,(-i)} - \hat{\mu}_t\right)^2,}
and the same form applies to the linked standard deviation
\eqn{\hat{\sigma}_t}.
The sandwich estimators use the item-additive estimating equations
of the closed-form methods (Robitzsch, 2023, 2024): with bread
\eqn{\mathbf{A}} and meat
\eqn{\mathbf{B}} the covariance of the linking parameters is
\eqn{\mathbf{V} = \mathbf{A}^{-1} \mathbf{B}
(\mathbf{A}^{-1})^{\top}}, and the linking error is its delta-method
projection onto the trend parameters. The bias-corrected variants
remove the part of the dispersion that is explained by the person
sampling of the item parameters and decompose the result into
standard error, linking error, and total error. All errors are
reported in the metric of the link's reference group.
}
\section{Which linking method supports which estimator}{
The refit jackknife and its bias-corrected variant are available for
every linking method and approach. The approximate jackknife (ajk)
covers the chain methods mgm, mm, haebara, and sl (sl with
\code{pow = 2}) and the pairwise Haberman variant (PHL). The
sandwich estimators and all remaining \code{"_bc"} variants cover
the chain methods mgm and mm and the pairwise Haberman variant.
With 3PL guessing parameters the approximate jackknife for haebara
and sl is unavailable (its closed forms are derived without
guessing). Use the refit jackknife there. For polytomous links the
refit jackknife and the approximate jackknife (chain mgm and mm,
PHL) are available. The sandwich and bias-corrected estimators are
derived for dichotomous models.
}
\examples{
# Identified parameters with a known trend and one DIF item:
mu_t <- c(0, 0.3, 0.6)
sigma_t <- c(1, 1.1, 1.2)
b_ref <- seq(-1.5, 1.5, length.out = 10)
ip <- data.frame(
  group = rep(1:3, each = 10),
  item = rep(sprintf("I\%02d", 1:10), 3),
  a = rep(sigma_t, each = 10),
  b = c(b_ref, (b_ref - 0.3) / 1.1, (b_ref - 0.6) / 1.2)
)
ip$b[ip$group == 2 & ip$item == "I01"] <-
  ip$b[ip$group == 2 & ip$item == "I01"] + 0.5
lk <- link_groups(ip, approach = "chain", method = "mgm")
linking_error(lk, method = "jackknife")$le

# Bias-corrected sandwich linking error needs an item parameter
# covariance. A small diagonal matrix stands in for a real one here
# (in practice use calibrate(..., keep_vcov = TRUE)$vcov_ipars).
V <- diag(0.001, 2 * nrow(lk$ipars))
linking_error(lk, method = "sandwich_osw_bc", vcov = V)$le

}
\references{
Monseur, C., & Berezner, A. (2007). The computation of equating
errors in international surveys in education. Journal of Applied
Measurement, 8(3), 323-335.

Robitzsch, A. (2023). Linking error in the 2PL model. J, 6(1),
58-84. \doi{10.3390/j6010005}

Robitzsch, A. (2024). Estimation of standard error, linking error,
and total error for robust and nonrobust linking methods in the
two-parameter logistic model. Stats, 7(3), 592-612.
\doi{10.3390/stats7030036}

Robitzsch, A. (2025a). Linking error estimation in Haberman linking.
AppliedMath, 5(1), 7. \doi{10.3390/appliedmath5010007}

Robitzsch, A. (2025b). Linking error estimation in Stocking-Lord
linking. Foundations, 5(1), 2. \doi{10.3390/foundations5010002}
}
