# Review of lane adhess: quadrature = TRUE and importance fits, whose
# outer objective may have no random effects in TMB's sense.
#   REV_ARM=lane|base Rscript dev/adhess-rev-quad.R
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-rev-common.R")
set.seed(31)
n <- 300
xs <- as.numeric(Sys.getenv("REV_XS", "100"))
d <- data.frame(x = rnorm(n, 5 * xs, xs), g = factor(rep(1:30, length.out = n)))
d$y <- rpois(n, exp(0.3 + 0.4 * (d$x - 5 * xs) / xs + rnorm(30, 0, 0.5)[d$g]))
for (what in c("quadrature", "importance")) {
  r <- cap(if (what == "quadrature") {
    frm(y ~ x + (1 | g), family = poisson(), data = d, quadrature = TRUE)
  } else {
    frm(y ~ x + (1 | g), family = poisson(), data = d, importance = 200L)
  })
  fit <- r$v
  if (inherits(fit, "rev_err")) {
    cat(what, "ERROR", fit, "\n")
    next
  }
  at <- the_point(fit)
  V <- suppressWarnings(vcov(fit, full = TRUE))
  Hr <- rich_hess(fit$obj, at)
  e <- cov_err(V, solve(Hr))
  hc <- fit$cache$hessian_fixed
  h <- tryCatch(ns$fit_outer_hessian(fit), error = function(e) NULL)
  cat(sprintf(paste0("%-10s re %d | exact %s | vcov vs Richardson of obj$gr ",
                     "se %.2e cor %.2e | lost %s | warn %d\n"),
              what, length(fit$obj$env$random), format(isTRUE(h$exact)),
              e[["se"]], e[["cor"]],
              paste(names(sdr_of(fit)$se_lost %||% character(0)),
                    collapse = ","), length(r$w)))
  cat("  SE:", format(sqrt(diag(V)), digits = 8), "\n")
  if (length(r$w)) cat("  warn:", substr(r$w, 1, 120), sep = "\n  ")
}
