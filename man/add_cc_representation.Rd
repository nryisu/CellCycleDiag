\name{add_cc_representation}
\alias{add_cc_representation}
\title{Add Seurat scores and tricycle cyclic position}
\description{Add Seurat scores and tricycle cyclic position.}
\usage{
add_cc_representation(object, species, gname.type = "SYMBOL",
  assay = "RNA", s.features = NULL, g2m.features = NULL,
  seed = 1, overwrite = FALSE)
}
\arguments{
\item{object}{A Seurat object with normalized RNA data. For diagnose_cc_imbalance, a metadata data.frame is also accepted.}
\item{species}{Explicit species: "human" or "mouse".}
\item{gname.type}{Gene identifiers: "SYMBOL" or "ENSEMBL"; must match assay row names.}
\item{assay}{Assay containing normalized expression; defaults to "RNA".}
\item{s.features, g2m.features}{Non-overlapping phase gene sets matching assay row names. Default human SYMBOL sets come from Seurat cc.genes.updated.2019. Both must be supplied for mouse or ENSEMBL inputs; no automatic orthology/case conversion is performed.}
\item{seed}{Non-negative integer seed. The caller's random-number state is restored.}
\item{overwrite}{Whether to replace existing representation metadata (default FALSE).}
}
\details{Requires Seurat, SeuratObject, SingleCellExperiment, SummarizedExperiment and tricycle. Supply a previously normalized assay with one joined data layer. The function does not normalize or renormalize expression. Seurat CellCycleScoring computes scores and tricycle projects normalized data using its built-in reference. For ENSEMBL inputs, tricycle may require species annotation packages. Inspect representation quality biologically, particularly in weakly cycling data. At least two matched genes per phase are required; this is an input check, not a quality threshold.}
\value{A Seurat object with S.Score, G2M.Score, Phase, theta, cc_sin and cc_cos metadata. Provenance is stored in Misc(object, slot = "CellCycleDiag"). Existing identities and default assay are preserved.}
