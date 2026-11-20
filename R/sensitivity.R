.check_design <- function(meta, group.by, covariates, label) {
  d <- data.frame(group = factor(meta[[group.by]]))
  for (i in seq_along(covariates)) {
    v <- meta[[covariates[i]]]
    if (anyNA(v) || (is.numeric(v) && any(!is.finite(v))))
      stop("Missing/non-finite covariate: ", covariates[i], call. = FALSE)
    if (!is.numeric(v) && !is.factor(v) && !is.character(v) && !is.logical(v))
      stop("Unsupported covariate: ", covariates[i], call. = FALSE)
    if (length(unique(v)) < 2L)
      stop("Constant covariate: ", covariates[i], call. = FALSE)
    d[[paste0("v", i)]] <- if (is.character(v) || is.factor(v)) factor(v) else v
  }
  mm <- stats::model.matrix(~ ., data = d)
  if (qr(mm)$rank < ncol(mm))
    stop(label, " design is rank deficient: check confounding/collinear covariates.", call. = FALSE)
  if (nrow(mm) <= ncol(mm))
    stop(label, " design has no residual degrees of freedom.", call. = FALSE)
  invisible(TRUE)
}

.sig_genes <- function(deg, padj_cutoff, logfc_cutoff, min.pct) {
  # An empty result is not accepted as proof of no DE: FindMarkers may skip testing.
  if (!is.data.frame(deg) || !nrow(deg))
    stop("DE returned no tested genes; cannot estimate correction susceptibility.", call. = FALSE)
  required <- c("avg_log2FC", "p_val_adj", "pct.1", "pct.2")
  if (!all(required %in% names(deg)))
    stop("DE results must contain ", paste(required, collapse = ", "), call. = FALSE)
  if (is.null(rownames(deg)) || anyDuplicated(rownames(deg)))
    stop("DE results require unique gene row names.", call. = FALSE)
  .numeric_columns(deg, required)
  if (any(deg$p_val_adj < 0 | deg$p_val_adj > 1) ||
      any(deg$pct.1 < 0 | deg$pct.1 > 1 | deg$pct.2 < 0 | deg$pct.2 > 1))
    stop("DE probability/detection fraction outside [0, 1].", call. = FALSE)
  rownames(deg)[deg$p_val_adj < padj_cutoff &
    abs(deg$avg_log2FC) > logfc_cutoff & pmax(deg$pct.1, deg$pct.2) >= min.pct]
}

compare_deg_sets <- function(direct, adjusted) {
  a <- unique(direct)
  b <- unique(adjusted)
  n_union <- length(union(a, b))
  j <- if (n_union) length(intersect(a, b))/n_union else NA_real_
  data.frame(direct_n = length(a), adjusted_n = length(b),
    overlap_n = length(intersect(a, b)), removed_n = length(setdiff(a, b)),
    gained_n = length(setdiff(b, a)),
    removed_fraction = if (length(a)) length(setdiff(a, b))/length(a) else NA_real_,
    jaccard_similarity = j, jaccard_distance = 1-j)
}

.fit_three <- function(object, group.by, ident.1, ident.2, assay, models,
                       min.pct, features, seed, verbose, finder) {
  ans <- lapply(names(models), function(nm) {
    tryCatch(.with_seed(seed, finder(object = object, group.by = group.by,
      ident.1 = ident.1, ident.2 = ident.2, assay = assay,
      test.use = "MAST", latent.vars = models[[nm]], slot = "data",
      min.pct = min.pct, logfc.threshold = 0, features = features,
      base = 2, only.pos = FALSE, random.seed = seed, verbose = verbose)),
      error = function(e) stop(nm, " model failed: ", conditionMessage(e), call. = FALSE))
  })
  stats::setNames(ans, names(models))
}

run_cc_sensitivity <- function(object, group.by, ident.1, ident.2,
                                donor.var = NULL, covariates = NULL,
                                assay = "RNA", s_col = "S.Score",
                                g2m_col = "G2M.Score", theta_col = "theta",
                                n_bins = 12, origin = 0, min.pct = 0.1,
                                padj_cutoff = 0.05, logfc_cutoff = 0.25,
                                features = NULL, seed = 1, verbose = FALSE) {
  .need("Seurat")
  .need("MAST")
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object.", call. = FALSE)
  .number(min.pct, "min.pct", 0, 1)
  .number(padj_cutoff, "padj_cutoff", 0, 1)
  .number(logfc_cutoff, "logfc_cutoff", 0)
  .number(seed, "seed", 0, .Machine$integer.max, TRUE)
  .scalar_string(ident.1, "ident.1")
  .scalar_string(ident.2, "ident.2")
  if (!is.null(donor.var)) .scalar_string(donor.var, "donor.var")
  if (!is.null(covariates) && (!is.character(covariates) || anyNA(covariates) || any(!nzchar(covariates))))
    stop("covariates must be character column names.", call. = FALSE)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  x <- .select_groups(.metadata(object), group.by, c(ident.1, ident.2))
  if (any(x$n < 3L)) stop("At least three cells per group are required for DE.", call. = FALSE)
  diagnostic <- diagnose_cc_imbalance(object, group.by, ident.1, ident.2,
                                       s_col, g2m_col, theta_col, n_bins, origin)
  baseline <- unique(c(donor.var, covariates))
  if (length(setdiff(baseline, names(x$meta)))) stop("Missing baseline covariates.", call. = FALSE)
  if (any(baseline %in% c(group.by, s_col, g2m_col, theta_col, "cc_sin", "cc_cos")))
    stop("Baseline covariates must exclude group and cell-cycle representation columns.", call. = FALSE)
  # Restrict the object itself so all downstream models use exactly this cohort.
  object <- object[, rownames(x$meta)]
  if (!is.null(donor.var)) {
    x$meta[[donor.var]] <- factor(x$meta[[donor.var]])
    object[[donor.var]] <- x$meta[[donor.var]]
  }
  m <- .normalized_data(object, assay)
  if (!is.null(features)) {
    if (!is.character(features) || !length(features) || anyNA(features) ||
        anyDuplicated(features) || !all(features %in% rownames(m)))
      stop("features must be unique gene names present in the selected assay.", call. = FALSE)
  }
  # Recompute cyclic covariates from the selected theta, ignoring any stale columns.
  cyclic <- c(".CellCycleDiag_sin", ".CellCycleDiag_cos")
  if (any(cyclic %in% c(group.by, baseline, s_col, g2m_col, theta_col)))
    stop("Reserved internal covariate names are in use.", call. = FALSE)
  theta <- x$meta[[theta_col]] %% (2*pi)
  x$meta[[cyclic[1]]] <- sin(theta)
  x$meta[[cyclic[2]]] <- cos(theta)
  object <- SeuratObject::AddMetaData(object, x$meta[, cyclic, drop = FALSE])
  models <- list(direct = baseline, score_adjusted = c(baseline, s_col, g2m_col),
                 cyclic_adjusted = c(baseline, cyclic))
  for (nm in names(models)) .check_design(x$meta, group.by, models[[nm]], nm)
  tables <- .fit_three(object, group.by, ident.1, ident.2, assay, models,
                       min.pct, features, seed, verbose, Seurat::FindMarkers)
  genes <- lapply(tables, .sig_genes, padj_cutoff = padj_cutoff,
                  logfc_cutoff = logfc_cutoff, min.pct = min.pct)
  tested <- lapply(tables, rownames)
  if (!all(vapply(tested[-1], setequal, logical(1), tested[[1]])))
    stop("Models returned different tested-gene universes; inspect failed/skipped tests.", call. = FALSE)
  impact <- do.call(rbind, lapply(names(genes)[-1], function(nm)
    cbind(model = nm, compare_deg_sets(genes$direct, genes[[nm]]))))
  rownames(impact) <- NULL
  structure(list(diagnostic = diagnostic, deg_tables = tables, deg_sets = genes,
    susceptibility = impact, models = models,
    settings = list(assay = assay, test.use = "MAST", donor.var = donor.var,
      covariates = covariates, min.pct = min.pct, padj_cutoff = padj_cutoff,
      logfc_cutoff = logfc_cutoff, features = features, seed = seed,
      tested_genes = tested[[1]], cells = colnames(object)), call = match.call(),
    versions = .versions(c("Seurat", "SeuratObject", "MAST", "CellCycleDiag")),
    session_info = utils::sessionInfo()), class = "cc_sensitivity")
}
