\name{irtlink-plots}
\alias{irtlink-plots}
\alias{plot.irtlink_dif}
\alias{plot.irtlink_link}
\alias{plot.irtlink_trend}
\title{Plot methods for irtlink objects}
\usage{
\method{plot}{irtlink_dif}(x, ...)

\method{plot}{irtlink_link}(x, pars = c("b", "a"), ...)

\method{plot}{irtlink_trend}(x, ...)
}
\arguments{
\item{x}{An \code{irtlink_dif}, \code{irtlink_link}, or
\code{irtlink_trend} object.}

\item{pars}{Linking plot only. Plot difficulties (\code{"b"},
default) or discriminations (\code{"a"}).}

\item{...}{Passed on to \code{\link[graphics:plot.default]{graphics::plot()}}.}
}
\value{
\code{x}, invisibly.
}
\description{
One plot method per result type. For a differential item
functioning (DIF) detection result the
items are sorted by their statistic (maximum root mean square
deviation (RMSD) or likelihood-ratio test (LRT) chi-square)
and the detection threshold is drawn as a reference line, with
flagged items highlighted. For a linking result each pair of
successive groups gets one panel that plots the common items'
parameters in the first group against the parameters of the second
group transformed onto the first group's metric. Items without
DIF or item parameter drift (IPD) scatter around the identity
line and outliers stick out.
For a trend result the estimated group means are drawn against the
group index.
}
