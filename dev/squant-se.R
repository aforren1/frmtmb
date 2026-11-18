# Lane squant, item 4(c): Hessian and sandwich standard errors of
# quantile_elf() against the Monte Carlo spread of the estimate, and
# against quantreg::rq(), on three designs.
#
# Run: Rscript dev/squant-se.R <design> <nrep> <n>
#   design: homo (gaussian errors), hetero (sd 0.5 + x), skew (lognormal
#   errors); seeds 20261116 + rep. Writes dev/squant-out/se-<design>-<n>.tsv
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(quantreg)
})
args <- commandArgs(TRUE)
design <- args[1]
nrep <- as.integer(args[2])
n <- as.integer(args[3])
taus <- c(0.5, 0.9)
hs <- c(0.1, 0.25, 0.5, 1)
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")

gen <- function(n) {
  x <- runif(n)
  e <- switch(design,
              homo = rnorm(n),
              hetero = (0.5 + x) * rnorm(n),
              skew = rlnorm(n))
  data.frame(x = x, y = 1 + 2 * x + e, row = seq_len(n))
}
truth <- function(tau) {
  switch(design,
         homo = c(1 + qnorm(tau), 2),
         hetero = c(1 + 0.5 * qnorm(tau), 2 + qnorm(tau)),
         skew = c(1 + qlnorm(tau), 2))
}

out <- file.path("dev/squant-out", sprintf("se-%s-%d.tsv", design, n))
cat("design\tn\trep\ttau\tmethod\th\tconv\tb0\tb1\tse_hess0\tse_hess1\t",
    "se_sand0\tse_sand1\ttrue0\ttrue1\twarnings\n", sep = "", file = out)
for (rep in seq_len(nrep)) {
  set.seed(20261116 + rep)
  dd <- gen(n)
  for (tau in taus) {
    tr <- truth(tau)
    r <- rq(y ~ x, tau = tau, data = dd)
    se <- tryCatch(summary(r, se = "nid")$coefficients[, 2],
                   error = function(e) c(NA, NA))
    row <- function(...) {
      cat(paste(c(...), collapse = "\t"), "\n", sep = "", file = out,
          append = TRUE)
    }
    row(design, n, rep, tau, "rq", NA, 0, coef(r), se, NA, NA, tr, "")
    for (h in hs) {
      warns <- character(0)
      res <- tryCatch(withCallingHandlers({
        f <- frm(bf(y ~ x, quantile = tau),
                 family = quantile_elf(smoothing = h), data = dd)
        b <- fixef(f)[, "Estimate"]
        sh <- sqrt(diag(vcov(f)))[1:2]
        ss <- sqrt(diag(vcov_cluster(f, ~ row, "CR0")))[1:2]
        c(f$opt$convergence, b, sh, ss)
      }, warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) rep(NA, 7))
      row(design, n, rep, tau, "elf", h, res, tr,
          gsub("[\t\n]", " ", paste(substr(warns, 1, 80), collapse = " | ")))
    }
  }
}
cat("done", design, n, nrep, "\n")
