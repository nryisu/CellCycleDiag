# Local validation — 2026-09-24

Environment: Windows x86_64, R 4.2.2. Library contains base/recommended packages.

## Completed

- `R CMD build --no-build-vignettes CellCycleDiag_package`: PASS.
- Windows binary built with `R CMD INSTALL --build`, installed into a separate
  workspace library via `install.packages(..., type="win.binary")`, and loaded
  successfully. The six exports and core suite passed against that installed binary.
- `R CMD check --no-manual --ignore-vignettes --no-build-vignettes CellCycleDiag_0.1.0.tar.gz`
  with `_R_CHECK_FORCE_SUGGESTS_=false`: **0 errors, 0 warnings, 1 NOTE**.
- Check passed installation, namespace/load/unload, syntax/static checks, Rd
  documentation and formal-argument consistency, examples, and base-R tests.
- Base tests cover known KS/TV values; empirical KS agreement with stats::ks.test;
  ties, symmetry, identical/separated distributions; theta wrapping, boundaries,
  origin; invalid/missing values; group selection; result classes; plotting;
  DEG cutoffs and empty-set semantics; rank-deficient designs; RNG restoration;
  and three-model dispatch with a test double, including failure propagation.
- Real-data compatibility: 50 mESC and 50 hESC v2 saved sampling records were
  joined to their corresponding saved metadata. Both package metrics were compared
  with the original local scripts on each of these 100 comparisons.
  Maximum absolute difference: KS components/D_score = 1.110223e-16;
  D_theta = 0. No differential-expression model was refitted in this audit.

## Not verified locally

- The one check NOTE lists unavailable optional dependencies: Seurat,
  SeuratObject, SingleCellExperiment, SummarizedExperiment, tricycle, MAST,
  testthat, knitr and rmarkdown. Dependency index retrieval was unavailable.
- Real Seurat/tricycle representation and Seurat/MAST DE integration tests are
  supplied in tests/testthat/test-integration.R but could not execute locally.
  Passing dispatch tests is not a substitute for executing these backends.
- PDF manual generation and rebuilding the knitr vignette were skipped.
  A prebuilt offline HTML guide is included; it was generated from the README,
  not from execution of the full Seurat workflow.
- GitHub CI configuration is included but has not run. No public repository or
  CRAN/Bioconductor submission has been created.

See RELEASE_CHECKLIST.md for release prerequisites.
