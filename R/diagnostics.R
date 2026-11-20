# Exact empirical KS distance, including ties; no p-value or asymptotics.
.ks_distance <- function(a, b) {
  at <- sort(unique(c(a, b)))
  max(abs(stats::ecdf(a)(at) - stats::ecdf(b)(at)))
}

calc_D_score <- function(meta, group_col, s_col = "S.Score",
                         g2m_col = "G2M.Score", group_levels = NULL) {
  x <- .select_groups(meta, group_col, group_levels)
  if (identical(s_col, g2m_col)) stop("Score columns must be distinct.", call. = FALSE)
  .numeric_columns(x$meta, c(s_col, g2m_col))
  ds <- .ks_distance(x$meta[[s_col]][x$a], x$meta[[s_col]][!x$a])
  dg <- .ks_distance(x$meta[[g2m_col]][x$a], x$meta[[g2m_col]][!x$a])
  cbind(.comparison(x), D_score_S = ds, D_score_G2M = dg, D_score = (ds + dg)/2)
}

calc_D_theta <- function(meta, group_col, theta_col = "theta", n_bins = 12,
                         group_levels = NULL, origin = 0) {
  .number(n_bins, "n_bins", 2, 100000, TRUE)
  .number(origin, "origin")
  x <- .select_groups(meta, group_col, group_levels)
  .numeric_columns(x$meta, theta_col)
  # Keep the original manuscript implementation: right-closed cut bins.
  angle <- (x$meta[[theta_col]] - origin) %% (2*pi)
  bins <- cut(angle, seq(0, 2*pi, length.out = n_bins + 1L),
              include.lowest = TRUE, labels = FALSE, right = TRUE)
  pa <- tabulate(bins[x$a], nbins = n_bins)/x$n[1]
  pb <- tabulate(bins[!x$a], nbins = n_bins)/x$n[2]
  cbind(.comparison(x), D_theta = sum(abs(pa-pb))/2,
        n_theta_bins = n_bins, theta_origin = origin)
}

diagnose_cc_imbalance <- function(object, group.by, ident.1, ident.2,
                                  s_col = "S.Score", g2m_col = "G2M.Score",
                                  theta_col = "theta", n_bins = 12, origin = 0) {
  .scalar_string(ident.1, "ident.1")
  .scalar_string(ident.2, "ident.2")
  x <- .select_groups(.metadata(object), group.by, c(ident.1, ident.2))
  sc <- calc_D_score(x$meta, group.by, s_col, g2m_col, x$groups)
  th <- calc_D_theta(x$meta, group.by, theta_col, n_bins, x$groups, origin)
  structure(list(comparison = .comparison(x),
    metrics = cbind(sc[, c("D_score_S", "D_score_G2M", "D_score"), drop = FALSE],
                    th[, "D_theta", drop = FALSE]),
    settings = list(group.by = group.by, s_col = s_col, g2m_col = g2m_col,
                    theta_col = theta_col, n_bins = n_bins, origin = origin,
                    excluded_cells = x$excluded),
    call = match.call(), package_version = "0.1.0"), class = "cc_diagnostic")
}
