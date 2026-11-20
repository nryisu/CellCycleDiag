add_cc_representation <- function(object, species, gname.type = "SYMBOL",
                                  assay = "RNA", s.features = NULL,
                                  g2m.features = NULL, seed = 1,
                                  overwrite = FALSE) {
  .need("Seurat")
  .need("SingleCellExperiment")
  .need("SummarizedExperiment")
  .need("tricycle")
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object.", call. = FALSE)
  species <- match.arg(species, c("human", "mouse"))
  gname.type <- match.arg(gname.type, c("SYMBOL", "ENSEMBL"))
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite))
    stop("overwrite must be TRUE or FALSE.", call. = FALSE)
  target <- c("S.Score", "G2M.Score", "Phase", "theta", "cc_sin", "cc_cos")
  existing <- intersect(target, names(object[[]]))
  if (!overwrite && length(existing))
    stop("Existing representation columns: ", paste(existing, collapse = ", "),
         ". Set overwrite = TRUE explicitly to replace them.", call. = FALSE)
  m <- .normalized_data(object, assay)
  if (is.null(s.features) || is.null(g2m.features)) {
    if (species != "human" || gname.type != "SYMBOL")
      stop("Supply both phase gene sets in the assay's identifiers for mouse/ENSEMBL data.", call. = FALSE)
    e <- new.env(parent = emptyenv())
    utils::data("cc.genes.updated.2019", package = "Seurat", envir = e)
    if (is.null(s.features)) s.features <- e$cc.genes.updated.2019$s.genes
    if (is.null(g2m.features)) g2m.features <- e$cc.genes.updated.2019$g2m.genes
  }
  match_genes <- function(g, label) {
    if (!is.character(g) || anyNA(g) || !length(g) || any(!nzchar(g)))
      stop(label, " must be a non-empty character vector without missing values.", call. = FALSE)
    kept <- intersect(g, rownames(m))
    if (length(kept) < 2L) stop("Fewer than two matched genes for ", label, call. = FALSE)
    if (length(kept) < length(unique(g)))
      warning(label, ": matched ", length(kept), "/", length(unique(g)), " genes.", call. = FALSE)
    kept
  }
  sf <- match_genes(s.features, "s.features")
  gf <- match_genes(g2m.features, "g2m.features")
  if (length(intersect(sf, gf))) stop("S and G2M gene sets must not overlap.", call. = FALSE)
  original_assay <- SeuratObject::DefaultAssay(object)
  object <- .with_seed(seed, {
    scored <- Seurat::CellCycleScoring(object, s.features = sf, g2m.features = gf,
                                     set.ident = FALSE, assay = assay, seed = seed)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = m))
    sce <- tricycle::project_cycle_space(sce, species = species, gname.type = gname.type)
    sce <- tricycle::estimate_cycle_position(sce)
    theta <- SummarizedExperiment::colData(sce)$tricyclePosition
    names(theta) <- colnames(sce)
    theta <- theta[colnames(scored)]
    if (!is.numeric(theta) || any(!is.finite(theta)))
      stop("tricycle did not return a finite angle for every cell.", call. = FALSE)
    theta <- theta %% (2*pi)
    scored <- SeuratObject::AddMetaData(scored, data.frame(theta = theta,
      cc_sin = sin(theta), cc_cos = cos(theta), row.names = colnames(scored)))
    scored
  })
  SeuratObject::DefaultAssay(object) <- original_assay
  provenance <- list(species = species, gname.type = gname.type, assay = assay,
    s.features = sf, g2m.features = gf, seed = seed, call = match.call(),
    versions = .versions(c("Seurat", "SeuratObject", "tricycle")))
  SeuratObject::Misc(object, slot = "CellCycleDiag") <- provenance
  object
}
