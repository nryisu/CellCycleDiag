\name{calc_D_theta}
\alias{calc_D_theta}
\title{Binned cyclic-position imbalance}
\description{Binned cyclic-position imbalance.}
\usage{
calc_D_theta(meta, group_col, theta_col = "theta", n_bins = 12,
  group_levels = NULL, origin = 0)
}
\arguments{
\item{meta}{Metadata data.frame; selected cells must have finite numeric representations.}
\item{group_col}{Name of the metadata column defining the comparison.}
\item{theta_col}{Name of the numeric angular coordinate column, in radians.}
\item{n_bins}{Integer number of equal-width angular bins (default 12; minimum 2).}
\item{group_levels}{Two distinct group labels, or NULL to require exactly two observed, non-missing groups.}
\item{origin}{Common angular origin in radians, default zero.}
}
\details{Angles are shifted by origin and wrapped modulo 2*pi. D_theta is one half of the sum of absolute differences in bin proportions. For manuscript compatibility, bins are right-closed and include zero in the first bin; 2*pi wraps to zero. This binned distance depends on bin count and origin and is not a rotation-invariant circular distance. Missing or non-finite selected angles cause an error.}
\value{One-row data.frame with group labels, cell counts, D_theta, n_theta_bins and theta_origin.}
\examples{
meta <- data.frame(condition = rep(c("Control", "Treatment"), each = 4),
                   S.Score = c(0, 0, 1, 1, 0, 1, 1, 2),
                   G2M.Score = c(0, 0, 0, 0, 1, 1, 1, 1),
                   theta = c(0, 0, pi, pi, 0, pi, pi, pi))
calc_D_theta(meta, "condition")
}
