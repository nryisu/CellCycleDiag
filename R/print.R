print.cc_diagnostic <- function(x, ..., digits = 3) {
  cat("CellCycleDiag analysis\n\nComparison\n")
  z <- x$comparison
  cat(sprintf("  %s: %s cells\n  %s: %s cells\n", z$group_A,
    format(z$n_A, big.mark = ",", trim = TRUE), z$group_B,
    format(z$n_B, big.mark = ",", trim = TRUE)))
  cat("\nPre-correction cell-cycle-state imbalance\n")
  m <- x$metrics
  values <- c(m$D_score, m$D_score_S, m$D_score_G2M, m$D_theta)
  labels <- c("  D_score ", "    KS(S) ", "    KS(G2M)", "  D_theta ")
  cat(paste0(labels, "  ", formatC(values, digits = digits, format = "f"), collapse = "\n"), "\n")
  cat("\nThese values quantify between-group cell-cycle-state\n",
      "distribution imbalance and do not determine whether\n",
      "cell-cycle correction is biologically required.\n", sep = "")
  invisible(x)
}

print.cc_sensitivity <- function(x, ..., digits = 3) {
  print(x$diagnostic, digits = digits)
  cat("\nCorrection susceptibility (relative to Direct)\n")
  print(x$susceptibility, digits = digits, row.names = FALSE)
  cat("\nDEG-set changes measure analytical sensitivity, not biological accuracy.\n")
  invisible(x)
}
