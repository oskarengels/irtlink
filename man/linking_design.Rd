\name{linking_design}
\alias{linking_design}
\title{Inspect the linking design of a data set list}
\usage{
linking_design(data_list)
}
\arguments{
\item{data_list}{List of data.frames or matrices, one per group, with item
names as column names.}
}
\value{
An object of class \code{irtlink_design}: list with \code{matrix}
(logical, items x groups), \code{successive_common} (integer vector),
and \code{n_groups}.
}
\description{
Builds the item-by-group incidence matrix (which item is administered
in which group) and counts the common items between successive
groups. Useful for checking the item overlap before calibration.
}
\examples{
d <- sim_trend_data(n_groups = 3, N = 20, I = 10)
linking_design(d)
}
