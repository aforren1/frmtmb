root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
sub1 <- function(path, old, new) {
  f <- file.path(root, path)
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  n <- length(gregexpr(old, txt, fixed = TRUE)[[1]])
  if (!grepl(old, txt, fixed = TRUE)) stop("no match in ", path)
  if (n != 1L) stop(n, " matches in ", path)
  writeLines(strsplit(sub(old, new, txt, fixed = TRUE), "\n",
                      fixed = TRUE)[[1]], f, useBytes = TRUE)
  cat("edited", path, "\n")
}
p <- "tests/testthat/test-scale-contract.R"

# The first yardstick was wrong, and the run said so: the mean and the
# median of a lognormal differ by about ONE standard error of the
# prediction on this design (median gap 0.9, not the 5 asserted), so
# the prediction's own error cannot separate them. The quantity that
# does is the ratio itself, whose standard error follows from the
# standard error of log(sigma) that summary() reports.
sub1(p,
paste0("  # then the discrimination, against the run's OWN standard error: the\n",
       "  # gap between the mean and the median must be large compared with\n",
       "  # how well the fit knows the prediction at all\n",
       "  se <- predict(f, type = \"link\", se.fit = TRUE)$se.fit\n",
       "  gap <- abs(fitted(f) - exp(mu)) / (fitted(f) * se)\n",
       "  expect_gt(stats::median(gap), 5)\n"),
paste0("  # then the discrimination. The prediction's own standard error\n",
       "  # cannot do it: measured on this design the mean and the median\n",
       "  # sit 0.9 of one apart. The quantity that separates them is the\n",
       "  # RATIO, exp(sigma^2 / 2), whose standard error follows by the\n",
       "  # delta method from the standard error of log(sigma) that\n",
       "  # summary() reports.\n",
       "  ratio <- stats::median(fitted(f) / exp(mu))\n",
       "  se_ls <- summary(f)$coefficients$sigma[1, 2]\n",
       "  se_ratio <- exp(sg^2 / 2) * sg^2 * se_ls\n",
       "  expect_gt((ratio - 1) / se_ratio, 5)\n"))

sub1(p,
paste0("  se <- stats::coef(summary(gs)$coefficients$mu)\n",
       "  se_int <- summary(gs)$coefficients$mu[1, 2]\n"),
       "  se_int <- summary(gs)$coefficients$mu[1, 2]\n")
cat("DONE\n")
