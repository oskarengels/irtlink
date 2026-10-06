\name{link_groups}
\alias{link_groups}
\title{Linking across groups}
\usage{
link_groups(
  x,
  approach = c("chain", "joint", "joint_restricted"),
  method = NULL,
  variant = c("simultaneous", "pairwise"),
  ref = NULL,
  ...
)
}
\arguments{
\item{x}{An \code{irtlink_calib} object or a data.frame with columns
\verb{group, item, a, b} (identified item parameters), optionally with a
column \code{c} of 3PL guessing parameters or centered threshold
columns \code{tau1, tau2, ...} for polytomous items.}

\item{approach}{\code{"chain"} estimates the linear transformation
between successive groups and accumulates the steps. \code{"joint"}
estimates the transformations of all groups in one fit.
\code{"joint_restricted"} repeats the joint fit on the growing window
of the first t groups and takes each group's result from the fit up
to that group.}

\item{method}{\code{"mgm"} (mean-geometric mean), \code{"mm"}
(mean-mean), \code{"haberman"}, \code{"haebara"}, or \code{"sl"}
(Stocking-Lord). The default is \code{"mgm"} for the chain and
\code{"haberman"} for the joint approaches. \code{"mm"} and
\code{"mgm"} are available for the chain approach only.}

\item{variant}{Form of the joint methods. \code{"simultaneous"}
(default) fits one criterion over all groups at once.
\code{"pairwise"} fits the criterion over all group pairs with
common items. For \code{method = "haberman"} the pairwise variant is
pairwise Haberman linking (PHL, Robitzsch, 2025).}

\item{ref}{Reference group whose metric carries the reported trend
(\code{mu = 0}, \code{sigma = 1} for this group). Default \code{NULL}
uses the first group. The estimation is unchanged and the trend is
re-expressed in the reference metric. \code{\link[=linking_error]{linking_error()}}
reports all errors in that metric, so the reference group has
linking error 0 and the first group then carries its own.}

\item{...}{Method options. \code{pow} is the Lp loss power (0, 0.25,
0.5, 1, 2) of \code{"haberman"}, \code{"haebara"}, and \code{"sl"}.
\code{weights} selects the quadrature weight scheme of
\code{"haebara"} and \code{"sl"} (\code{"uniform"},
\code{"normal_0.5"}, \code{"normal_1"}, \code{"normal_2"}).
\code{use_intercepts} makes \code{"haberman"} link the item
intercepts instead of the difficulties. \code{anchor} restricts the
linking to a character vector of item names (default: all common
items). \code{type} and \code{theta_metric} select the alignment
direction and the quadrature metric of the chain \code{"haebara"}
and \code{"sl"} criteria. \code{item_weights} weights the items of
the pairwise Haberman variant. \code{eps} overrides the smoothing
constant of the Lp loss.}
}
\value{
An object of class \code{irtlink_link}: a list with \code{trend}
(data.frame \verb{group, mu, sigma}), \code{steps} (chain approach
only, one row per successive pair), \code{method}, \code{approach},
\code{ipars} (the normalized item parameters), \code{anchor},
\code{le} (filled by \code{\link[=linking_error]{linking_error()}}), and \code{call}.
}
\description{
Places separately calibrated groups on one common scale. The chain
approach accumulates one step per pair of successive groups
(e.g., Kolen & Brennan, 2014),
\deqn{\mu_{t+1} = \mu_t + \sigma_t \, \hat{\mu}_{step},\qquad
\sigma_{t+1} = \sigma_t \, \hat{\sigma}_{step}.}
The joint approaches estimate all transformations in one fit. The
methods are the closed-form mean-mean and mean-geometric mean
linking, Haberman linking (Haberman, 2009) with the robust Lq loss
(Robitzsch, 2020) and a pairwise variant (Robitzsch, 2025), and the
response function methods of Haebara (1980) and Stocking and Lord
(1983). The Haebara criterion minimizes the discrepancy of the item
response functions over the common items,
\deqn{\sum_{i} \int \rho\big(P_i(\theta) -
P_i^{*}(\theta)\big)\, w(\theta)\, d\theta,}
and the Stocking-Lord criterion the discrepancy of the test
response functions,
\deqn{\int \rho\Big(\sum_i P_i(\theta) - \sum_i
P_i^{*}(\theta)\Big)\, w(\theta)\, d\theta,}
where \eqn{P_i^{*}} is the second group's item response function
under the candidate transformation, \eqn{w} the quadrature weights,
and \eqn{\rho(x) = |x|^p} the Lp loss, smoothed as
\eqn{(x^2 + \varepsilon)^{p/2}}. With 3PL item parameters the
response functions are the full \eqn{c_i + (1 - c_i) \,
\mathrm{logit}^{-1}(a_i (\theta - b_i))} curves (Kolen & Brennan,
2014). The guessing parameter enters the criteria but is itself
unaffected by the transformation. The simultaneous joint Haebara
fit re-estimates free item parameters and requires \eqn{c_i = 0},
its pairwise variant supports the 3PL.

For polytomous item parameters (generalized partial credit model,
threshold columns \code{tau1, tau2, ...}) the chain Haebara
criterion compares the category response functions and the chain
Stocking-Lord criterion the expected score functions (Kim & Lee,
2006). The pairwise joint criteria use the same polytomous response
functions. Mean-mean, mean-geometric mean, Haberman, and pairwise
Haberman linking operate on the discriminations and item locations
and need no further input. The simultaneous joint Haebara criterion
does not yet support polytomous items.
}
\examples{
# Identified parameters on the group-1 metric with a known trend:
# mu = (0, 0.3, 0.6), sigma = (1, 1.1, 1.2)
mu_t <- c(0, 0.3, 0.6)
sigma_t <- c(1, 1.1, 1.2)
b_ref <- seq(-1.5, 1.5, length.out = 10)
ip <- data.frame(
  group = rep(1:3, each = 10),
  item = rep(sprintf("I\%02d", 1:10), 3),
  a = rep(sigma_t, each = 10),
  b = c(b_ref, (b_ref - 0.3) / 1.1, (b_ref - 0.6) / 1.2)
)
link_groups(ip, approach = "chain", method = "mgm")
link_groups(ip, approach = "joint", method = "haberman",
            variant = "pairwise")

}
\references{
Haebara, T. (1980). Equating logistic ability scales by a weighted
least squares method. Japanese Psychological Research, 22(3),
144-149. \doi{10.4992/psycholres1954.22.144}

Haberman, S. J. (2009). Linking parameter estimates derived from an
item response model through separate calibrations (Research Report
RR-09-40). ETS. \doi{10.1002/j.2333-8504.2009.tb02197.x}

Kim, S., & Lee, W.-C. (2006). An extension of four IRT linking
methods for mixed-format tests. Journal of Educational Measurement,
43(1), 53-76. \doi{10.1111/j.1745-3984.2006.00004.x}

Kolen, M. J., & Brennan, R. L. (2014). Test equating, scaling, and
linking: Methods and practices (3rd ed.). Springer.
\doi{10.1007/978-1-4939-0317-7}

Robitzsch, A. (2020). Lp loss functions in invariance alignment and
Haberman linking with few or many groups. Stats, 3(3), 246-283.
\doi{10.3390/stats3030019}

Robitzsch, A. (2025). Linking error estimation in Haberman linking.
AppliedMath, 5(1), 7. \doi{10.3390/appliedmath5010007}

Stocking, M. L., & Lord, F. M. (1983). Developing a common metric in
item response theory. Applied Psychological Measurement, 7(2),
201-210. \doi{10.1177/014662168300700208}
}
