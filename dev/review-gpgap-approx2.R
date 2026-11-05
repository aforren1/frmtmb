# Review of lane gpgap, item 7: approximate gp() on UNSORTED positions, brms's
# own program with no rewriting. Plain gp(x, k = ) has only the sorted row in
# test-brms-likelihood.R. Reads the check's numbers; seed 8101.
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
H <- new.env(parent = asNamespace("frmtmb"))
for (h in c("helper-brms.R", "helper-stan.R")) {
  sys.source(file.path("tests/testthat", h), envir = H)
}
set.seed(8101)
n <- 80
d <- data.frame(x = sample(seq(0, 10, by = 0.25), n, TRUE),
                z = runif(n, 0, 5),
                f = factor(sample(c("a", "b"), n, TRUE)))
d$y <- sin(d$x) + 0.5 * cos(d$z) + rnorm(n, 0, 0.3)
stopifnot(is.unsorted(d$x))
rows <- list(
  k10 = y ~ gp(x, k = 10),
  k10_by = y ~ gp(x, by = f, k = 8),
  k5_2d = y ~ gp(x, z, k = 5),
  k5_2d_aniso = y ~ gp(x, z, k = 5, iso = FALSE),
  k10_grF = y ~ gp(x, k = 10, gr = FALSE))
for (nm in names(rows)) {
  fo <- rows[[nm]]
  out <- NULL
  r <- tryCatch({
    fit <- suppressWarnings(frm(do.call(frmtmb::bf, list(fo)), family = gaussian(),
                                data = d))
    with_reporter("silent", out <- H$brms_lp_check(brms::bf(fo), gaussian(), d, fit,
                                                   joint = TRUE))
    sprintf("const %.3e grad %.2e tol %.1e conv %d", out$measured_const, out$max_grad,
            1e-6 * abs(out$ours), fit$opt$convergence)
  }, error = function(e) paste("ERROR", substr(conditionMessage(e), 1, 150)))
  cat(sprintf("APPROX unsorted %-14s %s\n", nm, r))
}
