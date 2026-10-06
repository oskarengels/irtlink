\docType{package}
\name{irtlink-package}
\alias{irtlink}
\alias{irtlink-package}
\title{irtlink: IRT Linking and Linking Error Estimation Across Time Points
and Groups}
\description{
Estimates trends in the mean and standard deviation of latent ability
distributions across the time points of educational assessments, and
more generally ability distribution differences across linked groups,
with unidimensional item response theory (IRT) models (1PL, 2PL, 3PL),
together with the linking error that differential item functioning
(DIF) induces in the estimates.
}
\details{
The main function groups are:
\itemize{
\item Calibration: \code{\link[=calibrate]{calibrate()}} fits the groups on the internal
marginal maximum likelihood EM engine (Bock & Aitkin, 1981;
Hanson, 1998). It
covers separate calibration, concurrent calibration with full or
partial invariance, fixed item parameter calibration, and
regularized calibration with a smooth Bayesian information
criterion (SBIC) following O'Neill and Burke (2023) and Robitzsch
(2024b), for the 1PL, 2PL, and 3PL, with optional person sampling
weights (\code{pweights}).
\code{\link[=person_scores]{person_scores()}} returns expected a posteriori (EAP) and
weighted likelihood (WLE, Warm, 1989) person estimates.
\item Linking: \code{\link[=link_groups]{link_groups()}} selects the approach (chain, joint,
or restricted joint), the method (mean-mean, mean-geometric mean,
Haberman, 2009, with robust Lp loss, Haebara, 1980, and
Stocking-Lord, 1983), the simultaneous or pairwise form of the
joint methods, and the reference group. Concurrent calibration
likewise comes as chain, joint, or restricted joint composition
(via \code{\link[=estimate_trend]{estimate_trend()}}).
\item Linking error: \code{\link[=linking_error]{linking_error()}} estimates the item sampling
error of a linked trend by the delete-one item jackknife, its
bias-corrected and approximate (one-step) variants, and expected,
observed, and bias-corrected sandwich estimators (Robitzsch, 2023,
2024a, 2025a, 2025b). A \code{cluster} argument turns the jackknife
units into whole testlets. \code{\link[=dependent_vcov]{dependent_vcov()}} provides the joint
item parameter covariance for dependent person samples.
\item DIF diagnostics: \code{\link[=detect_dif]{detect_dif()}} flags items via root mean square
deviation (RMSD) statistics with fixed or data-driven cutoffs (von
Davier & Bezirhan, 2023) or per-item likelihood-ratio tests (Thissen,
Steinberg, & Wainer, 1988), with optional iterative purification.
\item Convenience: \code{\link[=estimate_trend]{estimate_trend()}} runs the full
calibrate-detect-link-error pipeline, \code{\link[=sim_trend_data]{sim_trend_data()}} simulates
longitudinal response data (uniform and non-uniform DIF, dependent
person samples), and \code{\link[=linking_design]{linking_design()}} inspects the item overlap
between groups.
}
}
\section{Model support}{

The 1PL and 2PL run through the complete pipeline, with every
calibration strategy, every linking method, DIF detection, and every
linking error estimator. The 3PL is supported for separate
calibration on the internal EM engine, with Bayes-modal estimation
of the guessing parameter. The Haebara and Stocking-Lord criteria
then compare the full 3PL response functions (Kolen & Brennan,
2014), while mean-mean, mean-geometric mean, and Haberman linking
operate on the discrimination and difficulty parameters by
construction. Linking errors under guessing are available through
the delete-one item jackknife and its bias-corrected variant. The
closed-form approximate jackknife for Haebara and Stocking-Lord is
derived without guessing.

Polytomous items are supported through the generalized partial
credit model (GPCM, Muraki, 1992) and the partial credit model
(PCM, Masters, 1982), with separate calibration, chain and pairwise
joint linking by every method (Haebara over the category response
functions and Stocking-Lord over the expected score functions, Kim
& Lee, 2006), joint Haberman and pairwise Haberman linking, linking
errors by the refit jackknife and the approximate jackknife for
mean-mean and mean-geometric mean chains, and EAP person scores.
}

\section{Related packages}{

Linking and equating methods are also provided by \pkg{sirt} and
\pkg{TAM} (Robitzsch), \pkg{equateIRT} (Battauz), \pkg{kequate}
(Andersson, Br\enc{ä}{a}nberg, & Wiberg), and \pkg{plink} (Weeks).
General IRT estimation is provided by \pkg{mirt} (Chalmers).
irtlink differs in its focus on linking error estimation for trend
estimates.
}

\references{
Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood
estimation of item parameters: Application of an EM algorithm.
Psychometrika, 46(4), 443-459. \doi{10.1007/BF02293801}

Haebara, T. (1980). Equating logistic ability scales by a weighted
least squares method. Japanese Psychological Research, 22(3),
144-149. \doi{10.4992/psycholres1954.22.144}

Haberman, S. J. (2009). Linking parameter estimates derived from an
item response model through separate calibrations (Research Report
RR-09-40). ETS. \doi{10.1002/j.2333-8504.2009.tb02197.x}

Hanson, B. A. (1998). IRT parameter estimation using the EM
algorithm. Technical note.

Kim, S., & Lee, W.-C. (2006). An extension of four IRT linking
methods for mixed-format tests. Journal of Educational Measurement,
43(1), 53-76. \doi{10.1111/j.1745-3984.2006.00004.x}

Kolen, M. J., & Brennan, R. L. (2014). Test equating, scaling, and
linking: Methods and practices (3rd ed.). Springer.
\doi{10.1007/978-1-4939-0317-7}

Masters, G. N. (1982). A Rasch model for partial credit scoring.
Psychometrika, 47(2), 149-174. \doi{10.1007/BF02296272}

Muraki, E. (1992). A generalized partial credit model: Application
of an EM algorithm. Applied Psychological Measurement, 16(2),
159-176. \doi{10.1177/014662169201600206}

O'Neill, M., & Burke, K. (2023). Variable selection using a smooth
information criterion for distributional regression models.
Statistics and Computing, 33, 71. \doi{10.1007/s11222-023-10204-8}

Robitzsch, A. (2023). Linking error in the 2PL model. J, 6(1),
58-84. \doi{10.3390/j6010005}

Robitzsch, A. (2024a). Estimation of standard error, linking error,
and total error for robust and nonrobust linking methods in the
two-parameter logistic model. Stats, 7(3), 592-612.
\doi{10.3390/stats7030036}

Robitzsch, A. (2024b). Smooth information criterion for regularized
estimation of item response models. Algorithms, 17(4), 153.
\doi{10.3390/a17040153}

Robitzsch, A. (2025a). Linking error estimation in Haberman linking.
AppliedMath, 5(1), 7. \doi{10.3390/appliedmath5010007}

Robitzsch, A. (2025b). Linking error estimation in Stocking-Lord
linking. Foundations, 5(1), 2. \doi{10.3390/foundations5010002}

Stocking, M. L., & Lord, F. M. (1983). Developing a common metric in
item response theory. Applied Psychological Measurement, 7(2),
201-210. \doi{10.1177/014662168300700208}

Thissen, D., Steinberg, L., & Wainer, H. (1988). Use of item response
theory in the study of group difference in trace lines. In H. Wainer
& H. Braun (Eds.), Test validity (pp. 147-170). Lawrence Erlbaum.

von Davier, M., & Bezirhan, U. (2023). A robust method for detecting
item misfit in large-scale assessments. Educational and Psychological
Measurement, 83(4), 740-765. \doi{10.1177/00131644221105819}

Warm, T. A. (1989). Weighted likelihood estimation of ability in
item response theory. Psychometrika, 54(3), 427-450.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/oskarengels/irtlink}
  \item Report bugs at \url{https://github.com/oskarengels/irtlink/issues}
}

}
\author{
\strong{Maintainer}: Oskar Engels \email{engels@leibniz-ipn.de}

}
\keyword{internal}
