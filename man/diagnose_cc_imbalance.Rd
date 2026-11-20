\name{diagnose_cc_imbalance}
\alias{diagnose_cc_imbalance}
\title{Diagnose between-group cell-cycle-state imbalance}
\description{Diagnose between-group cell-cycle-state imbalance.}
\usage{
diagnose_cc_imbalance(object, group.by, ident.1, ident.2,
  s_col = "S.Score", g2m_col = "G2M.Score", theta_col = "theta",
  n_bins = 12, origin = 0)
}
\arguments{
\item{object}{A Seurat object with normalized RNA data. For diagnose_cc_imbalance, a metadata data.frame is also accepted.}
\item{group.by}{Name of the metadata column defining the comparison.}
\item{ident.1, ident.2}{Distinct character group labels. Other groups and missing group labels are excluded explicitly.}
\item{s_col, g2m_col}{Names of the two score columns.}
\item{theta_col}{Name of the numeric angular coordinate column, in radians.}
\item{n_bins}{Integer number of equal-width angular bins (default 12; minimum 2).}
\item{origin}{Common angular origin in radians, default zero.}
}
\details{Both representations must already exist. Use add_cc_representation to compute them, or supply existing metadata. The same selected cells are used for both metrics. No automatic missing-value deletion, phase inference, correction decision, or risk classification is performed.}
\value{An S3 cc_diagnostic list with comparison, metrics, settings (including excluded_cells), call and package_version. The input object is unchanged.}
\examples{
meta <- data.frame(condition = rep(c("Control", "Treatment"), each = 4),
                   S.Score = c(0, 0, 1, 1, 0, 1, 1, 2),
                   G2M.Score = c(0, 0, 0, 0, 1, 1, 1, 1),
                   theta = c(0, 0, pi, pi, 0, pi, pi, pi))
res <- diagnose_cc_imbalance(meta, "condition", "Control", "Treatment")
res
plot_cc_diagnostic(res)
}
