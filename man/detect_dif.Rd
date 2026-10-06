\name{detect_dif}
\alias{detect_dif}
\title{Detect differential item functioning}
\usage{
detect_dif(
  calib,
  method = c("rmsd", "lrt"),
  rmsd_method = c("parameter", "data", "both"),
  cutoff_type = c("fixed", "data_driven"),
  cutoff = 0.05,
  tau = 2.7,
  flag_on = c("parameter", "data"),
  link_method = "mgm",
  procedure = c("one_step", "iterative_forward"),
  max_iter = 7,
  min_anchor = 3,
  alpha = 0.05,
  alpha_adjust = c("none", "bonferroni"),
  dif_type = c("both", "uniform", "nonuniform"),
  ...
)
}
\arguments{
\item{calib}{An \code{irtlink_calib} object (separate calibration,
\code{keep_models = TRUE}).}

\item{method}{\code{"rmsd"} (default) or \code{"lrt"} (likelihood-ratio test).}

\item{rmsd_method}{\code{"parameter"} (default), \code{"data"}, or \code{"both"} (RMSD).}

\item{cutoff_type}{\code{"fixed"} (default) or \code{"data_driven"} (RMSD).}

\item{cutoff}{Fixed RMSD cutoff(s). A vector is allowed with
\code{procedure = "one_step"}. Default \code{0.05}.}

\item{tau}{Robust z cutoff for \code{cutoff_type = "data_driven"}. Default \code{2.7}.}

\item{flag_on}{Which estimator drives flagging when \code{rmsd_method = "both"}.}

\item{link_method}{Internal linking method for the parameter-based RMSD.}

\item{procedure}{\code{"one_step"} (default) or \code{"iterative_forward"}. The LRT is
always one-step.}

\item{max_iter}{Maximum purification iterations. Default \code{7}.}

\item{min_anchor}{Minimum anchor items to retain. Default \code{3}.}

\item{alpha}{Significance level for \code{method = "lrt"}. Default \code{0.05}.}

\item{alpha_adjust}{Multiple-comparison adjustment for the LRT p-values:
\code{"none"} (default) or \code{"bonferroni"}.}

\item{dif_type}{For \code{method = "lrt"}, the kind of DIF/IPD tested
per item. \code{"both"} (default) frees discrimination and
difficulty with df = 2 per extra group, \code{"uniform"} the
difficulty only (df = 1), and \code{"nonuniform"} the discrimination
only (df = 1). Ignored for \code{method = "rmsd"}.}

\item{...}{Passed on to the internal fitting functions.}
}
\value{
An object of class \code{irtlink_dif}: a list with \code{flagged} (names
of the flagged items), \code{anchor} (retained items), \code{rmsd} (the driver RMSD table,
\code{NULL} for the LRT), \code{rmsd_parameter}/\code{rmsd_data} (per-estimator tables,
possibly \code{NULL}), \code{lrt} (the per-item LRT table \verb{item, chisq, df, p, p_adj} for \code{method = "lrt"}, else \code{NULL}), \code{flagged_by_cutoff},
\code{linkable}, \code{history} (per-iteration, or \code{NULL}), \code{converged}, and the
call settings (\code{method}, \code{dif_type}, \code{alpha}, etc.).
}
\description{
Flags items showing differential item functioning (DIF) across
groups, known in longitudinal assessments as item parameter drift
(IPD). Two detection methods are available. The root mean square
deviation (RMSD) method (\code{method = "rmsd"}) measures the
deviation between an item's group-specific response function and
its invariance-constrained reference function,
\deqn{RMSD_i = \sqrt{\int (P_{i}(\theta) - P_{i}^{ref}(\theta))^2
\, w(\theta) \, d\theta},}
with either a parameter-based (population) estimator from linked
separate calibrations or a data-based (empirical) estimator computed
from a full-invariance multiple group fit (see K\enc{ö}{oe}hler,
Robitzsch, & Hartig, 2020, and Tijmstra et al., 2020, for the two
estimators). The data-based tables also report bias-corrected RMSD
values that remove the positive finite sample bias of the sample
RMSD (K\enc{ö}{oe}hler et al., 2020), using the analytic estimator
T4 of Robitzsch (2025). Items are flagged with fixed cutoffs or with the
data-driven robust z cutoff of von Davier and Bezirhan (2023),
based on the median absolute deviation (MAD),
\deqn{z_i = \frac{RMSD_i - \mathrm{median}(RMSD)}{1.4826 \cdot
\mathrm{MAD}(RMSD)} > \tau,}
in one step or with an iterative forward purification in which
flagged items leave the anchor and stay excluded. Operational
large-scale assessments flag items at RMSD cutoffs around 0.12
(e.g., OECD, 2024), and cutoff choice is discussed in
K\enc{ö}{oe}hler et al. (2020) and von Davier and Bezirhan (2023).
In trend estimation simulations, the fixed cutoff 0.05 performed
well across conditions, 0.03 lost some precision in small samples,
0.08 was too lenient under strong balanced drift, and for the
data-driven cutoff \code{tau = 2.7} kept the better efficiency
while \code{tau = 1.7} with iterative purification minimized bias
in extreme conditions (Engels, L\enc{ü}{ue}dtke, & Robitzsch,
2026).
The likelihood-ratio test (\code{method = "lrt"}, Thissen,
Steinberg, & Wainer, 1988) splits each item across groups and
compares the deviance of a full-invariance multiple group model to
that of the item-split model, with all other items as anchor (Wang
& Yeh, 2003). \code{dif_type} chooses whether to test uniform
(difficulty), non-uniform (discrimination), or both kinds of
DIF/IPD.
}
\details{
Both methods require a separate calibration with \code{keep_models = TRUE}
(the per-group responses are reconstructed from the stored models).
At least \code{min_anchor} items always remain in the anchor (see
Kopf, Zeileis, & Strobl, 2015, for a framework of anchor selection
methods).
}
\examples{
set.seed(1)
d <- sim_trend_data(n_groups = 2, N = 400, I = 12, dif = "unbalanced",
                    dif_pct = 0.2, dif_effect = 0.8)
cal <- calibrate(d)
dif <- detect_dif(cal, rmsd_method = "parameter", cutoff = 0.1)
dif$flagged
dif$anchor

}
\references{
Engels, O., L\enc{ü}{ue}dtke, O., & Robitzsch, A. (2026). Estimating
trends with differential item functioning: A comparison of five
IRT-based approaches. Educational and Psychological Measurement.
\doi{10.1177/00131644251408818}

K\enc{ö}{oe}hler, C., Robitzsch, A., & Hartig, J. (2020). A
bias-corrected RMSD item fit statistic: An evaluation and comparison
to alternatives. Journal of Educational and Behavioral Statistics,
45(3), 251-273. \doi{10.3102/1076998619890566}

Kopf, J., Zeileis, A., & Strobl, C. (2015). A framework for anchor
methods and an iterative forward approach for DIF detection. Applied
Psychological Measurement, 39(2), 83-103.
\doi{10.1177/0146621614544195}

OECD (2024). PISA 2022 technical report. OECD Publishing.
\doi{10.1787/01820d6d-en}

Robitzsch, A. (2025). Bias-corrected root mean square deviation
estimators. Foundations, 5(4), 36. \doi{10.3390/foundations5040036}

Thissen, D., Steinberg, L., & Wainer, H. (1988). Use of item response
theory in the study of group difference in trace lines. In H. Wainer
& H. Braun (Eds.), Test validity (pp. 147-170). Lawrence Erlbaum.

Tijmstra, J., Bolsinova, M., Liaw, Y.-L., Rutkowski, L., & Rutkowski,
D. (2020). Sensitivity of the RMSD for detecting item-level misfit in
low-performing countries. Journal of Educational Measurement, 57(4),
566-583. \doi{10.1111/jedm.12263}

von Davier, M., & Bezirhan, U. (2023). A robust method for detecting
item misfit in large-scale assessments. Educational and Psychological
Measurement, 83(4), 740-765. \doi{10.1177/00131644221105819}

Wang, W.-C., & Yeh, Y.-L. (2003). Effects of anchor item methods on
differential item functioning detection with the likelihood ratio
test. Applied Psychological Measurement, 27(6), 479-498.
\doi{10.1177/0146621603259902}
}
