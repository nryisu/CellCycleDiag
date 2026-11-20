\name{plot_cc_diagnostic}
\alias{plot_cc_diagnostic}
\title{Plot descriptive imbalance or correction susceptibility}
\description{Plot descriptive imbalance or correction susceptibility.}
\usage{
plot_cc_diagnostic(x, type = c("imbalance", "susceptibility"),
  main = NULL, col = c("#2166AC", "#B35806"))
}
\arguments{
\item{x}{A cc_diagnostic or cc_sensitivity result.}
\item{type}{"imbalance", or "susceptibility" for a cc_sensitivity result.}
\item{main}{Optional plot title.}
\item{col}{Bar colours.}
}
\details{No risk bands, recommended thresholds or yes/no correction decisions are displayed. Undefined empty-set comparisons are labelled explicitly.}
\value{Invisibly returns the plotted data.frame; draws a base-R bar plot on the current graphics device.}
\examples{
meta <- data.frame(condition = rep(c("Control", "Treatment"), each = 4),
                   S.Score = c(0, 0, 1, 1, 0, 1, 1, 2),
                   G2M.Score = c(0, 0, 0, 0, 1, 1, 1, 1),
                   theta = c(0, 0, pi, pi, 0, pi, pi, pi))
res <- diagnose_cc_imbalance(meta, "condition", "Control", "Treatment")
res
plot_cc_diagnostic(res)
}
