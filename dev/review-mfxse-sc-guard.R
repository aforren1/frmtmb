# Review of lane mfxse, short check: m1's guard and cost claim, the row check
# under a user vcov =, and a function hypothesis on a fit.
# Data seed 7411 (mre_data(); 500 levels x 40 rows for the cost).
#   Rscript dev/review-mfxse-sc-guard.R
.libPaths(c("C:/Users/adf44/source/r/wt-mfxse-lib", "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({library(frmtmb); library(marginaleffects)})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
quiet <- function(e) suppressWarnings(suppressMessages(e))
ws <- function(e) {
  out <- character(0)
  v <- withCallingHandlers(tryCatch(e, error = function(err) paste("ERROR", conditionMessage(err))),
    warning = function(w) {
      out <<- c(out, paste0(class(w)[1], ": ", substr(conditionMessage(w), 1, 110)))
      invokeRestart("muffleWarning")
    })
  list(v = v, w = out)
}
set.seed(7411)
J <- 12L; m <- 20L; n <- J * m
d <- data.frame(x = rnorm(n), f = factor(sample(c("a", "b"), n, TRUE)),
                g = factor(rep(seq_len(J), each = m)))
d$cnt <- rpois(n, exp(0.5 + 0.4 * d$x + 0.3 * (d$f == "b") + rnorm(J, 0, 0.8)[d$g]))
f1 <- quiet(frm(cnt ~ x + f + (1 | g), family = poisson(), data = d))
nd <- d[c(1, 50, 120), ]
cat("\n== the guard (options(frmtmb.latres_dense_rows = 5), 12 inner rows)\n")
p0 <- quiet(predictions(f1, newdata = nd, vcov = ~ g))
old <- options(frmtmb.latres_dense_rows = 5L)
r <- ws(predictions(f1, newdata = nd, vcov = ~ g))
options(old)
held <- quiet(predictions(f1, newdata = nd, vcov = get_vcov(f1, vcov = ~ g)[1:3, 1:3]))
cat("carried SE ", sprintf("%.4f", p0$std.error), "\n")
cat("guarded SE ", if (is.character(r$v)) r$v else sprintf("%.4f", r$v$std.error), "\n")
cat("coef-only  ", sprintf("%.4f", held$std.error), "\n")
cat("warnings:", paste(r$w, collapse = "\n          "), "\n")
old <- options(frmtmb.latres_dense_rows = 5L)
r2 <- ws(predictions(f1, newdata = nd))
options(old)
cat("default vcov under the same option:", sprintf("%.4f", r2$v$std.error), " warnings", length(r2$w), "\n")
cat("\n== the row check under a user vcov = (fresh fit, nothing armed)\n")
f2 <- quiet(frm(cnt ~ x + f + (1 | g), family = poisson(), data = d))
ndn <- data.frame(x = c(0, 1, 2), f = "a", g = factor(c("3", "new", "new2")))
for (lab in c("default", "~ g", "vcov(full = TRUE)")) {
  r <- ws(switch(lab,
    default = predictions(f2, newdata = ndn, allow_new_levels = TRUE),
    "~ g" = predictions(f2, newdata = ndn, allow_new_levels = TRUE, vcov = ~ g),
    "vcov(full = TRUE)" = predictions(f2, newdata = ndn, allow_new_levels = TRUE,
                                      vcov = vcov(f2, full = TRUE))))
  k <- sum(grepl("^frmtmb_interop_se_part", r$w))
  cat(sprintf("  %-18s frmtmb_interop_se_part warnings %d  SE %s\n", lab, k,
              paste(sprintf("%.3f", r$v$std.error), collapse = " ")))
}
cat("\n== a function hypothesis on a fit (seed 7411, (1 | g))\n")
ap <- avg_predictions(f2)
hf <- quiet(hypotheses(f2, hypothesis = function(x) data.frame(term = "avg", estimate = avg_predictions(x)$estimate)))
old <- quiet(avg_predictions(f2, vcov = get_vcov(f2)[1:3, 1:3]))
cat(sprintf("  avg_predictions SE %.4f; hypotheses(fit, function) SE %.4f; coefficient-only (0.79.1 convention) %.4f\n",
            ap$std.error, hf$std.error, old$std.error))
cat("\n== m1's cost claim, 500 levels\n")
J <- 500L; m <- 40L; n <- J * m
set.seed(7411)
db <- data.frame(x = rnorm(n), g = factor(rep(seq_len(J), each = m)))
db$y <- rpois(n, exp(0.5 + 0.4 * db$x + rnorm(J, 0, 0.8)[db$g]))
fb <- quiet(frm(y ~ x + (1 + x | g), family = poisson(), data = db))
el <- function(e) { t0 <- proc.time()[["elapsed"]]; force(e); proc.time()[["elapsed"]] - t0 }
Q <- ns$joint_precision(fb); rr <- fb$obj$env$random
cat(sprintf("  inner rows %d\n", length(rr)))
cat(sprintf("  Matrix::solve(Qrr) dense       %.2f s\n", el(as.matrix(Matrix::solve(Q[rr, rr])))))
cat(sprintf("  vcov_cluster(~ g, full = TRUE) %.2f s\n", el(quiet(vcov_cluster(fb, ~ g, full = TRUE)))))
cat(sprintf("  get_vcov(vcov = ~ g)           %.2f s\n", el(quiet(get_vcov(fb, vcov = ~ g)))))
cat(sprintf("  predictions(2 rows, vcov = ~ g) %.2f s\n", el(quiet(predictions(fb, newdata = db[c(1, n), ], vcov = ~ g)))))
