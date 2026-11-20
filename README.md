# CellCycleDiag

An R package for **between-group cell-cycle-state imbalance** and optional
**correction susceptibility** analysis. 

Cell-cycle variation, between-group imbalance and correction susceptibility
are distinct concepts. No high/medium/low risk labels, automatic correction
thresholds or claims that cyclic adjustment is superior are provided.

## Install locally

```r
install.packages("CellCycleDiag_0.1.0.tar.gz", repos = NULL, type = "source")
library(CellCycleDiag)
```

The diagnostic and plotting functions use only base/recommended R packages.
For the full Seurat workflow install optional dependencies in your R environment:

```r
install.packages(c("Seurat", "BiocManager"))
BiocManager::install(c("tricycle", "MAST", "SingleCellExperiment",
                       "SummarizedExperiment"))
```

Use a Bioconductor release compatible with your R version. ENSEMBL projection
may additionally need `org.Hs.eg.db` or `org.Mm.eg.db`.

### From GitHub (recommended)

```r
install.packages("remotes")
remotes::install_github("nryisu/CellCycleDiag", force = TRUE, build_vignettes = TRUE)
library(CellCycleDiag)
```

## Runnable metadata example

```r
meta <- data.frame(condition = rep(c("Control", "Treatment"), each = 4),
                   S.Score = c(0, 0, 1, 1, 0, 1, 1, 2),
                   G2M.Score = c(0, 0, 0, 0, 1, 1, 1, 1),
                   theta = c(0, 0, pi, pi, 0, pi, pi, pi))
res <- diagnose_cc_imbalance(meta, "condition", "Control", "Treatment")
res
plot_cc_diagnostic(res)
```

## Seurat workflow

```r
# seu: your quality-controlled, log-normalized Seurat object.
# NormalizeData and, for split v5 assays, JoinLayers should be done beforehand.
seu <- add_cc_representation(seu, species = "human", assay = "RNA")
res <- diagnose_cc_imbalance(seu, group.by = "condition",
                             ident.1 = "Control", ident.2 = "Treatment")
res
plot_cc_diagnostic(res)

sens <- run_cc_sensitivity(seu, group.by = "condition",
                          ident.1 = "Control", ident.2 = "Treatment",
                          donor.var = "donor")
sens
sens$susceptibility
sens$deg_tables$score_adjusted
plot_cc_diagnostic(sens, type = "susceptibility")
saveRDS(sens, "cc_sensitivity.rds")
```

For mouse or ENSEMBL identifiers, explicitly supply validated `s.features`
and `g2m.features` matching assay row names. Existing scores/theta can be used
directly without `add_cc_representation()`. It refuses overwriting by default.
Use `overwrite = TRUE` only to deliberately recompute those columns.

## Six exported functions

| Function | Purpose |
|---|---|
| `add_cc_representation()` | Seurat S/G2M scores and tricycle theta |
| `calc_D_score()` | Mean of two empirical KS distances |
| `calc_D_theta()` | Total variation of binned theta distributions |
| `diagnose_cc_imbalance()` | Unified `cc_diagnostic` result |
| `run_cc_sensitivity()` | Direct, score-adjusted and cyclic-adjusted MAST DE |
| `plot_cc_diagnostic()` | Descriptive imbalance/susceptibility plots |

`D_score = (KS(S) + KS(G2M))/2`.
`D_theta = sum(abs(p_bin - q_bin))/2`, with 12 equal-width bins and origin 0.
Angles wrap modulo 2*pi. Bins are right-closed (include zero in the first bin),
matching the original manuscript code. Bin count/origin affect this approximation.

All selected cells must have finite scores and angles; there is no silent
pairwise deletion. Counts and excluded groups are recorded. Diagnostics do not
return a significance test and do not determine whether correction is needed.

DE adjustment is performed through MAST covariates, not by regressing the
expression matrix before DE. All models include the same baseline covariates
(including donor when supplied), use the same cells, min.pct filter and genes.
The default DEG rule is adjusted P < 0.05, absolute log2FC > 0.25, and maximum
detection fraction >= 0.1. For the manuscript's gradient analysis set
`logfc_cutoff = 0.10`. Seurat supplies Bonferroni-adjusted P values and marginal
fold changes; fold changes are not adjusted effect estimates.

Jaccard distance is 1 - intersection/union. Both-empty sets yield NA;
one-empty sets yield 1. Removed fraction is NA when Direct has no DEGs.
Failures stop the analysis instead of being converted to zero DEGs.
Rank deficiency (including complete donor/condition confounding) is rejected.
Donor fixed effects do not replace biological replication or implement pseudobulk.

## Provenance and migration

This package refactors the local CellCycleDiag manuscript scripts. The original
scripts and manuscript are preserved. Low-level metric interfaces retain
`meta`, `group_col`, and `group_levels`; DE and plotting use the new result types.
No original simulation, phase classifier or PCA/UMAP helper is exported.
Results retain calls, settings and, for sensitivity analysis, dependency versions
and session information. Representation provenance is in
`SeuratObject::Misc(seu, slot = "CellCycleDiag")`.

## Development and publication

Run `R CMD build CellCycleDiag_package` and `R CMD check CellCycleDiag_0.1.0.tar.gz`.
Full vignette/check tooling requires `testthat`, `knitr` and `rmarkdown`.
Base-R tests run even without testthat. Integration tests run automatically when
their optional dependencies are installed. See VALIDATION.md for actual local coverage.

The DESCRIPTION maintainer email is an explicit placeholder (`maintainer@example.org`)
and must be replaced before publication; it is not a verified author identity.
No GitHub repository has been published. After pushing this directory's contents
to a repository, `remotes::install_github("YOUR_ACCOUNT/CellCycleDiag")` can be used
with the actual account. See RELEASE_CHECKLIST.md.

## Method references

- [Seurat CellCycleScoring](https://satijalab.org/seurat/reference/cellcyclescoring)
- [Seurat FindMarkers](https://satijalab.org/seurat/reference/findmarkers)
- [tricycle vignette](https://bioconductor.org/packages/release/bioc/vignettes/tricycle/inst/doc/tricycle.html)
- Zheng et al. (2022), tricycle. https://doi.org/10.1186/s13059-021-02581-y
- Finak et al. (2015), MAST. https://doi.org/10.1186/s13059-015-0844-5
