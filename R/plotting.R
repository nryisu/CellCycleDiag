plot_cc_diagnostic <- function(x, type = c("imbalance", "susceptibility"),
                                main = NULL, col = c("#2166AC", "#B35806")) {
  type <- match.arg(type)
  if (!inherits(x, "cc_diagnostic") && !inherits(x, "cc_sensitivity"))
    stop("x must be a cc_diagnostic or cc_sensitivity result.", call. = FALSE)
  if (type == "imbalance") {
    d <- if (inherits(x, "cc_sensitivity")) x$diagnostic else x
    values <- c(D_score = d$metrics$D_score, D_theta = d$metrics$D_theta)
    if (is.null(main)) main <- "Between-group cell-cycle-state imbalance"
    out <- data.frame(metric = names(values), value = unname(values))
  } else {
    if (!inherits(x, "cc_sensitivity"))
      stop("Susceptibility plots require run_cc_sensitivity() output.", call. = FALSE)
    values <- stats::setNames(x$susceptibility$jaccard_distance,
                             c("Score-adjusted", "Cyclic-adjusted"))
    if (is.null(main)) main <- "Correction susceptibility: Jaccard distance"
    out <- data.frame(metric = names(values), value = unname(values))
  }
  # NA empty-set comparisons remain visibly undefined, rather than plotted as zero.
  positions <- graphics::barplot(replace(values, is.na(values), 0), ylim = c(0, 1.12),
    col = col, main = main, ylab = "Distance", border = NA)
  graphics::text(positions, ifelse(is.na(values), 0.04, values + 0.04),
    labels = ifelse(is.na(values), "undefined", formatC(values, digits = 3, format = "f")))
  invisible(out)
}
