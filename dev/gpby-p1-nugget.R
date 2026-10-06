# Punch round 1, m2, nugget decision: frmtmb's exact gp() with its
# nugget at 1e-6 of the correlation (the current value) against 1e-12
# (brms's jitter). For each, on the same data: does the fit converge and
# agree with the other, how ill-conditioned K is, whether K's Cholesky
# factor (prediction, kriging) still exists, and how rough a draw is
# between two positions 1e-9 apart. The value is the user's to choose;
# this measures, it does not change it.
# Usage: Rscript dev/gpby-p1-nugget.R <nugget>
nug <- as.numeric(commandArgs(TRUE)[1])
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
utils::assignInNamespace("gp_nugget", nug, "frmtmb")
cat("lib:", find.package("frmtmb"), "| gp_nugget", ns$gp_nugget, "\n")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-gpby"
src <- parse(file.path(wt, "tests/testthat/test-gp-by.R"))
for (e in src) if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e)
set.seed(23)
g2 <- local({
  grid <- expand.grid(x1 = seq(0, 4, by = 0.5), x2 = seq(0, 4, by = 0.5))
  pu <- as.matrix(grid[sample(nrow(grid), 50), ])
  idx <- sample(50, 80, replace = TRUE)
  data.frame(y = 1 + sin(pu[idx, 1]) * cos(pu[idx, 2]) +
               stats::rnorm(80, 0, 0.4), x1 = pu[idx, 1], x2 = pu[idx, 2])
})
set.seed(7)
dense <- data.frame(x = round(stats::runif(200, 0, 6), 2))
dense$y <- sin(dense$x) + stats::rnorm(200, 0, 0.3)
set.seed(8)
nd9 <- data.frame(x = sort(stats::runif(40, 0, 6)))
nd9 <- rbind(nd9, data.frame(x = nd9$x[1:10] + 1e-9))
nd9$y <- sin(nd9$x) + stats::rnorm(50, 0, 0.3)
cases <- list(
  list("gp(x), 60 rows", bf(y ~ gp(x)), gpby_data()),
  list("gp(x, by = f)", bf(y ~ gp(x, by = f)), gpby_data()),
  list("gp(x1, x2), 2-D", bf(y ~ gp(x1, x2)), g2),
  list("gp(x), 200 rows 0.01 apart", bf(y ~ gp(x)), dense),
  list("gp(x), 10 pairs 1e-9 apart", bf(y ~ gp(x)), nd9)
)
for (cs in cases) {
  lab <- cs[[1]]
  t <- system.time(fit <- tryCatch(suppressWarnings(
    frm(cs[[2]], data = cs[[3]])), error = function(e) e))[["elapsed"]]
  if (inherits(fit, "error")) {
    cat(sprintf("FIT %-28s | ERROR %s\n", lab, conditionMessage(fit)))
    next
  }
  th <- fit$estimates$theta
  bks <- Filter(function(b) b$covstruct == "gp", fit$frame$re_blocks)
  kap <- vapply(bks, function(bk) {
    thc <- th[bk$theta_idx]
    thc[1] <- 0
    kappa(ns$covstruct_registry$gp$vcov(thc, bk), exact = TRUE)
  }, 0)
  okchol <- vapply(bks, function(bk) {
    thc <- th[bk$theta_idx]
    thc[1] <- 0
    !inherits(tryCatch(chol(ns$covstruct_registry$gp$vcov(thc, bk)),
                       error = function(e) e), "error")
  }, NA)
  dg <- diagnose(fit, quiet = TRUE)
  cat(sprintf(paste0("FIT %-28s | logLik %.6f | conv %d (%s) | pdHess %s ",
                     "| max |grad| %.1e | max kappa(K) %.2e | chol(K) %s ",
                     "| %.1f s\n"), lab,
              as.numeric(logLik(fit)), as.integer(dg$convergence),
              dg$message, dg$pdHess, dg$max_grad, max(kap),
              paste(okchol, collapse = ","), t))
  cat("    theta", paste(sprintf("%.6f", th), collapse = " "), "\n")
  if (identical(lab, "gp(x), 60 rows")) {
    # two rows 1e-9 apart past the data and inside it: the difference of
    # one draw of the field between them, against the field's sd there
    nd <- data.frame(x = c(7.4, 7.4 + 1e-9, 2.55, 2.55 + 1e-9))
    lp <- fit$frame$linpreds[["y.mu"]]
    ed <- ns$lp_eta_design(fit, lp, nd, TRUE, FALSE)
    krig <- Filter(Negate(is.null), lapply(ed$sm_parts, `[[`, "krig"))[[1]]
    S <- ns$gp_krig_cov(krig)
    sdd <- c(sqrt(S[1, 1] + S[2, 2] - 2 * S[1, 2]),
             sqrt(S[3, 3] + S[4, 4] - 2 * S[3, 4]))
    cat(sprintf(paste0("    rows 1e-9 apart: sd of the draw's difference ",
                       "%.2e (x = 7.4), %.2e (x = 2.55); field sd there ",
                       "%.2e, %.2e; chol(K) in prediction %s\n"), sdd[1],
                sdd[2], sqrt(S[1, 1]), sqrt(S[3, 3]),
                if (is.null(krig$P)) "FAILED (solve fallback)" else "ok"))
    pr <- frm_linpred(fit, newdata = data.frame(x = c(2.55, 6.5, 7.5)),
                      se.fit = TRUE)
    cat("    se.fit at 2.55, 6.5, 7.5:",
        paste(sprintf("%.6f", pr$se.fit), collapse = " "), "\n")
  }
}
