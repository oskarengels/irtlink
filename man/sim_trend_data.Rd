\name{sim_trend_data}
\alias{sim_trend_data}
\title{Simulate longitudinal item response data}
\usage{
sim_trend_data(
  n_groups,
  N,
  I = 20,
  mu = NULL,
  sigma = NULL,
  overlap = 0.4,
  design = c("all", "successive"),
  dif = c("none", "balanced", "unbalanced", "random"),
  dif_pct = 0.3,
  dif_effect = 0.5,
  dif_sd = 0.3,
  dif_par = c("b", "a", "both"),
  dif_effect_a = 0.25,
  dif_sd_a = 0.15,
  person_overlap = 0,
  person_cor = 0.7,
  model = c("2PL", "1PL", "GPCM"),
  n_cat = 3
)
}
\arguments{
\item{n_groups}{Number of groups (at least 2).}

\item{N}{Persons per group.}

\item{I}{Items administered per group.}

\item{mu, sigma}{Numeric vectors of length \code{n_groups} with the latent means
\eqn{\mu_t} and standard deviations \eqn{\sigma_t} of the groups, in
group order. Defaults: \code{mu = 0.3 * (0:(n_groups-1))},
\code{sigma = rep(1, n_groups)}.}

\item{overlap}{Proportion of \code{I} used as common (link) items,
giving \code{round(I * overlap)} items.}

\item{design}{\code{"all"} shares one common item set across all
groups, \code{"successive"} places common items only between
neighboring groups.}

\item{dif}{Type of differential item functioning (DIF, in
longitudinal assessments item parameter drift): \code{"none"}
(default), \code{"balanced"}
(affected items shift in alternating directions), \code{"unbalanced"} (all
affected items shift by \code{+dif_effect}), or \code{"random"} (per-group N(0,
dif_sd) shifts). DIF is applied to link-relevant items
at groups >= 2. Note that \code{"random"} mode draws its shifts before
the responses, so it consumes the RNG stream differently from the
other modes.}

\item{dif_pct}{Proportion of link-relevant items affected by DIF. The
count is \code{round(n_link * dif_pct)} (round-half-to-even). Affected items
are chosen in alphabetical order, so for a successive design the DIF
concentrates in the earliest link blocks.}

\item{dif_effect}{Difficulty shift size for balanced/unbalanced DIF.}

\item{dif_sd}{Standard deviation of the per-group shifts for random
DIF.}

\item{dif_par}{Parameter carrying the DIF: \code{"b"} (default, uniform
DIF, difficulty shifts), \code{"a"} (non-uniform DIF, shifts on the log
discriminations, i.e. multiplicative), or \code{"both"}. \code{"a"}/\code{"both"}
require \code{model = "2PL"} and use \code{dif_effect_a}/\code{dif_sd_a} for the
log-discrimination shifts.}

\item{dif_effect_a, dif_sd_a}{Analogues of \code{dif_effect}/\code{dif_sd} on the
log-discrimination scale, used when \code{dif_par} is \code{"a"} or \code{"both"}.}

\item{person_overlap}{Proportion of the \code{N} persons per group that are
shared with the previous group (dependent person samples, default 0,
independent samples). Shared persons keep their identifier across
groups (row names and the \code{person_ids} attribute), so the group fits
can be passed to \code{\link[=dependent_vcov]{dependent_vcov()}}.}

\item{person_cor}{Latent correlation between the standardized
abilities of a shared person in consecutive groups (default 0.7,
only used when \code{person_overlap > 0}).}

\item{model}{\code{"2PL"} (default), \code{"1PL"} (all discriminations
fixed at 1), or \code{"GPCM"} (polytomous responses from the
generalized partial credit model).}

\item{n_cat}{Number of response categories per item for
\code{model = "GPCM"} (default 3). The step difficulties are
\code{b + tau_v} with centered, equally spaced \code{tau_v} between
-0.6 and 0.6, identical across items.}
}
\value{
List of \code{n_groups} data.frames (persons x administered items) with
attributes \code{true_mu}, \code{true_sigma}, \code{true_ipars} (data.frame
\verb{item, a, b}), \code{design}, \code{dif}, \code{dif_items} (names of the items with
DIF),
\code{dif_effects} (a pool-by-group matrix of injected difficulty shifts),
\code{dif_effects_a} (the analogous log-discrimination shifts, all zero
for \code{dif_par = "b"}), and \code{person_ids} (list of per-group person
identifiers, equal to the row names).
}
\description{
Generates dichotomous (2PL/1PL) or polytomous (GPCM) response data
for several groups, such
as time points, with a linking design. Group 1 is the reference
population: the first elements \code{mu[1]} and \code{sigma[1]} of the \code{mu} and
\code{sigma} vectors are its latent mean and standard deviation (typically
0 and 1). Trends are induced through the latent distribution.
}
\note{
Subsetting the returned list (e.g., \code{d[1:2]}) drops the attributes.
Keep the full object when the true values are needed.

In \code{random} mode the DIF shifts are drawn before the response
data, so the same seed yields different responses than
\code{dif = "none"}. The modes \code{balanced} and \code{unbalanced}
leave the response RNG stream unchanged.

The default item pool of ten discrimination and difficulty pairs,
recycled to the requested test length, follows the pool used in the
linking error simulation studies of Robitzsch (2023).
}
\references{
Robitzsch, A. (2023). Linking error in the 2PL model. J, 6(1),
58-84. \doi{10.3390/j6010005}
}
\examples{
d <- sim_trend_data(n_groups = 3, N = 100)
attr(d, "true_mu")

}
