\name{run_cc_sensitivity}
\alias{run_cc_sensitivity}
\title{Compare DE sensitivity to cell-cycle adjustment}
\description{Compare DE sensitivity to cell-cycle adjustment.}
\usage{
run_cc_sensitivity(object, group.by, ident.1, ident.2,
  donor.var = NULL, covariates = NULL, assay = "RNA",
  s_col = "S.Score", g2m_col = "G2M.Score", theta_col = "theta",
  n_bins = 12, origin = 0, min.pct = 0.1, padj_cutoff = 0.05,
  logfc_cutoff = 0.25, features = NULL, seed = 1, verbose = FALSE)
}
\arguments{
\item{object}{A Seurat object with normalized RNA data. For diagnose_cc_imbalance, a metadata data.frame is also accepted.}
\item{group.by}{Name of the metadata column defining the comparison.}
\item{ident.1, ident.2}{Distinct character group labels. Other groups and missing group labels are excluded explicitly.}
\item{donor.var}{Optional donor column, encoded categorically in all three models.}
\item{covariates}{Optional additional metadata covariates included in every model.}
\item{assay}{Assay containing normalized expression; defaults to "RNA".}
\item{s_col, g2m_col}{Names of the two score columns.}
\item{theta_col}{Name of the numeric angular coordinate column, in radians.}
\item{n_bins}{Integer number of equal-width angular bins (default 12; minimum 2).}
\item{origin}{Common angular origin in radians, default zero.}
\item{min.pct}{Minimum detection fraction in either group for testing and DEG selection (default 0.1).}
\item{padj_cutoff}{Strict upper cutoff for Seurat Bonferroni-adjusted P values (default 0.05).}
\item{logfc_cutoff}{Strict lower cutoff for absolute unadjusted log2 fold change (default 0.25).}
\item{features}{Optional unique assay gene names tested by all models, or NULL for all eligible features.}
\item{seed}{Non-negative integer seed. The caller's random-number state is restored.}
\item{verbose}{Whether Seurat prints progress.}
}
\details{Uses Seurat FindMarkers with MAST only, normalized data, base 2 and no logFC prefilter. Direct includes baseline covariates; score-adjusted adds S and G2M scores; cyclic-adjusted adds sin(theta) and cos(theta). All models use the same cells, feature eligibility and baseline covariates. Adjustment occurs in the DE model, not by residualizing the expression matrix. Reported log2 fold changes are marginal Seurat fold changes, not adjusted model coefficients. Donor is a fixed effect; this does not implement pseudobulk or remove cell-level pseudoreplication. Rank-deficient designs, missing covariates, failed models, invalid results and unequal tested-gene universes cause errors. Jaccard distance is 1 minus intersection/union. Both empty sets give NA; one empty set gives 1. Removed fraction is undefined (NA) if the Direct set is empty. DEG-set changes measure sensitivity, not correctness.}
\value{An S3 cc_sensitivity list containing diagnostic, deg_tables, deg_sets, susceptibility, models, settings, call, versions and session_info. The input object is unchanged.}
