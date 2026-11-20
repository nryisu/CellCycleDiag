\name{calc_D_score}
\alias{calc_D_score}
\title{Score-distribution imbalance}
\description{Score-distribution imbalance.}
\usage{
calc_D_score(meta, group_col, s_col = "S.Score",
  g2m_col = "G2M.Score", group_levels = NULL)
}
\arguments{
\item{meta}{Metadata data.frame; selected cells must have finite numeric representations.}
\item{group_col}{Name of the metadata column defining the comparison.}
\item{s_col, g2m_col}{Names of the two score columns.}
\item{group_levels}{Two distinct group labels, or NULL to require exactly two observed, non-missing groups.}
}
\details{D_score is the arithmetic mean of the two empirical two-sample Kolmogorov--Smirnov distances. Ties are evaluated exactly from empirical CDFs. No P value is computed. All distances lie in [0, 1]. Missing or non-finite selected scores cause an error.}
\value{One-row data.frame containing group labels, cell counts, D_score_S, D_score_G2M and D_score.}
\examples{
meta <- data.frame(condition = rep(c("Control", "Treatment"), each = 4),
                   S.Score = c(0, 0, 1, 1, 0, 1, 1, 2),
                   G2M.Score = c(0, 0, 0, 0, 1, 1, 1, 1),
                   theta = c(0, 0, pi, pi, 0, pi, pi, pi))
calc_D_score(meta, "condition")
}
