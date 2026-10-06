# irtlink
#### IRT Linking and Linking Error Estimation Across Time Points and Groups

irtlink uses unidimensional item response theory (IRT) models (1PL,
2PL, 3PL, and the partial credit and generalized partial credit
models for polytomous items) to estimate trends in the mean and
standard deviation of a latent ability distribution across time
points of an educational assessment. It also estimates ability
distribution differences across linked groups. Groups are connected
through common items; their person samples may be independent,
dependent, or drawn from different populations. The package
quantifies the linking error induced by differential item
functioning (DIF) in the estimates. For a tour of the package, see
`vignette("irtlink")`, and for the function overview, see
`?irtlink`.

If you encounter a bug or have a suggestion, please open a GitHub
issue or email me at engels@leibniz-ipn.de. The most useful reports
contain a small example that reproduces the problem, including a
minimal dataset (a `sim_trend_data()` call with a seed is usually
sufficient), the shortest script that triggers the issue on that
dataset, and the output of `sessionInfo()`.

#### GitHub version `irtlink` 0.1.0 (2026-10-06)

[![](https://img.shields.io/badge/github%20version-0.1.0-orange.svg)](https://github.com/oskarengels/irtlink)

The version hosted here can be installed from within R using:

```r
# install.packages("remotes")
remotes::install_github("oskarengels/irtlink")
```

or with `devtools`:

```r
# install.packages("devtools")
devtools::install_github("oskarengels/irtlink")
```
