.need <- function(package) {
  if (!requireNamespace(package, quietly = TRUE))
    stop("Install optional dependency '", package, "' to use this function.", call. = FALSE)
}

.scalar_string <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.character(x) || !nzchar(x))
    stop(name, " must be a non-empty character scalar.", call. = FALSE)
}

.number <- function(x, name, lower = -Inf, upper = Inf, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || x > upper || (integer && x != floor(x)))
    stop(name, " is outside its permitted range.", call. = FALSE)
}

.metadata <- function(object) {
  if (is.data.frame(object)) return(object)
  if (inherits(object, "Seurat")) {
    .need("SeuratObject")
    return(object[[]])
  }
  stop("object must be a Seurat object or metadata data.frame.", call. = FALSE)
}

.select_groups <- function(meta, group_col, group_levels = NULL) {
  if (!is.data.frame(meta) || anyDuplicated(names(meta)))
    stop("meta must be a data.frame with unique column names.", call. = FALSE)
  .scalar_string(group_col, "group_col")
  if (!group_col %in% names(meta)) stop("Missing group column: ", group_col, call. = FALSE)
  g <- as.character(meta[[group_col]])
  if (is.null(group_levels)) {
    if (anyNA(g)) stop("Missing group labels; specify two group_levels or clean labels.", call. = FALSE)
    group_levels <- unique(g)
  }
  group_levels <- as.character(group_levels)
  if (length(group_levels) != 2L || anyNA(group_levels) ||
      any(!nzchar(group_levels)) || anyDuplicated(group_levels))
    stop("Specify exactly two distinct, non-missing group labels.", call. = FALSE)
  keep <- !is.na(g) & g %in% group_levels
  ans <- meta[keep, , drop = FALSE]
  a <- g[keep] == group_levels[1L]
  if (!any(a) || all(a)) stop("Each selected group must contain cells.", call. = FALSE)
  list(meta = ans, a = a, groups = group_levels,
       n = c(sum(a), sum(!a)), excluded = sum(!keep))
}

.numeric_columns <- function(meta, columns) {
  for (nm in columns) {
    .scalar_string(nm, "column name")
    if (!nm %in% names(meta)) stop("Missing column: ", nm, call. = FALSE)
    if (!is.numeric(meta[[nm]]) || any(!is.finite(meta[[nm]])))
      stop(nm, " must contain finite numeric values for every selected cell.", call. = FALSE)
  }
}

.comparison <- function(x) data.frame(
  group_A = x$groups[1], group_B = x$groups[2],
  n_A = x$n[1], n_B = x$n[2], row.names = NULL
)

.normalized_data <- function(object, assay) {
  .need("SeuratObject")
  .scalar_string(assay, "assay")
  if (!assay %in% names(object)) stop("Assay not found: ", assay, call. = FALSE)
  a <- object[[assay]]
  if (inherits(a, "Assay5")) {
    layers <- SeuratObject::Layers(a)
    if (!"data" %in% layers || any(grepl("^data\\.", layers)))
      stop("A single joined normalized 'data' layer is required. Normalize/join layers first.", call. = FALSE)
    m <- SeuratObject::LayerData(a, layer = "data")
  } else {
    m <- SeuratObject::GetAssayData(object, assay = assay, slot = "data")
  }
  if (!nrow(m) || !ncol(m) || !setequal(colnames(m), colnames(object)))
    stop("Normalized assay data must cover every selected cell.", call. = FALSE)
  m[, colnames(object), drop = FALSE]
}

.with_seed <- function(seed, code) {
  .number(seed, "seed", 0, .Machine$integer.max, TRUE)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", old, envir = .GlobalEnv)
          else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  force(code)
}

.versions <- function(packages) {
  stats::setNames(vapply(packages, function(p)
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p))
    else NA_character_, character(1)), packages)
}
