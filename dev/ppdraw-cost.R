# What draw_edge_scan() costs, once per fit, against predict() itself.
# Usage: Rscript dev/ppdraw-cost.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
tm <- function(expr, reps) {
  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(reps)) force(expr())
  (proc.time()[["elapsed"]] - t0) / reps
}
set.seed(31)
J <- 100
n <- 3000
d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                g = factor(sample(seq_len(J), n, TRUE)),
                h = factor(sample(seq_len(30), n, TRUE)))
u <- matrix(stats::rnorm(2 * J, sd = c(0.5, 0.3)), J, 2, byrow = TRUE)
d$y <- stats::rpois(n, exp(0.3 + (0.4 + u[d$g, 2]) * d$x + 0.2 * d$z +
                             u[d$g, 1] + stats::rnorm(30, sd = 0.3)[d$h]))
fits <- list(
  glmm = frm(bf(y ~ x + z + (1 + x | g) + (1 | h)), family = poisson(),
             data = d),
  zinb = frm(bf(y ~ x + z + (1 | g), zi ~ z), family =
               zero_inflated_negbinomial(), data = d))
for (nm in names(fits)) {
  fit <- fits[[nm]]
  ds <- frmtmb:::fit_draw_space(fit)
  sd2 <- diag(ds$V)
  np <- length(fit$opt$par)
  ts <- tm(function() frmtmb:::draw_edge_scan(fit, ds$map, sd2), 5)
  tf <- tm(function() fit$obj$fn(fit$opt$par), 20)
  tp <- tm(function() {
    fit$cache$draw_edge_held <- integer(0)
    suppressMessages(predict(fit, ndraws = 200, summary = FALSE))
  }, 3)
  cat(sprintf(paste0("%s: np=%d, one objective %.4f s, scan %.3f s, ",
                     "predict(ndraws = 200) with the verdict cached ",
                     "%.3f s, scan / predict %.3f\n"),
              nm, np, tf, ts, tp, ts / tp))
}
