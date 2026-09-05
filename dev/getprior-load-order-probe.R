# Does get_prior() change its answer when frmtmb.sample is attached?
# Written 2026-09-04 as the demonstration behind the route= change; the
# regression test in tests/testthat/ is derived from it.
suppressMessages(library(frmtmb))
data(sleepstudy, package = "lme4")
f <- bf(Reaction ~ Days + (Days | Subject)) + gaussian()
fit <- frm(f, data = sleepstudy)
before_f <- capture.output(print(get_prior(f, data = sleepstudy)))
before_fit <- capture.output(print(get_prior(fit)))
cat("frmtmb.sample installed:", requireNamespace("frmtmb.sample", quietly = TRUE), "\n")
if (requireNamespace("frmtmb.sample", quietly = TRUE)) {
  suppressMessages(library(frmtmb.sample))
} else {
  # Walk up until extensions/frmtmb.sample is under foot, so the probe
  # runs from the repository root and from dev/ alike. The fixed
  # arithmetic this replaces assumed dev/ and resolved one level short
  # from the root; an installed frmtmb.sample skips the branch, which
  # is what kept that from being noticed.
  root <- normalizePath(".", winslash = "/")
  while (!dir.exists(file.path(root, "extensions", "frmtmb.sample")) &&
           !identical(dirname(root), root)) {
    root <- dirname(root)
  }
  suppressMessages(pkgload::load_all(
    file.path(root, "extensions", "frmtmb.sample"), quiet = TRUE))
}
after_f <- capture.output(print(get_prior(f, data = sleepstudy)))
after_fit <- capture.output(print(get_prior(fit)))
cat("\n== get_prior(formula) identical before/after attach:",
    identical(before_f, after_f), "==\n")
if (!identical(before_f, after_f)) {
  cat("--- before ---\n"); cat(before_f, sep = "\n")
  cat("--- after ---\n");  cat(after_f, sep = "\n")
}
cat("\n== get_prior(fit) identical before/after attach:",
    identical(before_fit, after_fit), "==\n")
if (!identical(before_fit, after_fit)) {
  cat("--- before ---\n"); cat(before_fit, sep = "\n")
  cat("--- after ---\n");  cat(after_fit, sep = "\n")
}
