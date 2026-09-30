# Reviewer, lane ordinal: post-fit methods on the new structures and
# links. fitted(scale = "linear") per threshold, prediction on new data,
# simulate, bootstrap refits, conditional_effects, emmeans, summary.
# Data seed 20261006. Output: dev/ordinal-rev-log-postfit.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(emmeans)})
set.seed(20261006)
n <- 300
d <- data.frame(x = rnorm(n, 1), z = rnorm(n),
                h = factor(sample(c("p", "q"), n, TRUE)))
u <- stats::rlogis(n) / exp(0.3 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -0.2) + (u > 0.8) + (u > 1.8) + (u > 2.8)
run <- function(lab, expr) {
  cat("\n##", lab, "\n")
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("  WARNING:", substr(conditionMessage(w), 1, 150), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))
  invisible(r)
}
fits <- list(
  equi_cs = frm(bf(y ~ x + cs(z), disc ~ 0 + z), data = d,
                family = sratio(threshold = "equidistant")),
  stz_cs = frm(y ~ x + cs(z), data = d, family = acat(threshold = "sum_to_zero")),
  cum_equi = frm(y ~ x, data = d, family = cumulative(threshold = "equidistant")),
  cum_stz = frm(y ~ x, data = d, family = cumulative(threshold = "sum_to_zero")),
  acat_probit = frm(y ~ x, data = d, family = acat("probit")),
  gr_equi = frm(y | thres(gr = h) ~ x, data = d,
                family = cratio(threshold = "equidistant"))
)
for (nm in names(fits)) {
  f <- fits[[nm]]
  cat("\n==========", nm, "\n")
  run("fitted(scale = 'linear') dim and rel. diff from X b (+ cs_k)", {
    fl <- fitted(f, scale = "linear")
    b <- fixef(f)
    eta <- d$x * b["x", "Estimate"]
    cat("  dim:", dim(fl), "\n")
    est <- if (length(dim(fl)) == 3L) fl[, "Estimate", ] else fl[, "Estimate"]
    if (is.matrix(est)) {
      cs <- b[grep("^z\\[", rownames(b)), "Estimate"]
      ref <- eta + outer(d$z, cs)
      cat("  per-threshold max rel diff:", max(abs(est - ref) / abs(ref)), "\n")
    } else {
      cat("  max rel diff:", max(abs(est - eta) / abs(eta)), "\n")
    }
  })
  run("fitted(newdata) == fitted rows", {
    nd <- d[c(3, 50, 200), ]
    a <- fitted(f, newdata = nd)[, "Estimate", ]
    b <- fitted(f)[c(3, 50, 200), "Estimate", ]
    cat("  max abs diff:", max(abs(a - b)), "\n")
  })
  run("simulate", { s <- simulate(f, nsim = 2, seed = 1)
    cat("  table:", table(unlist(s)), "\n") })
  run("predict", { p <- predict(f, ndraws = 50); cat("  dim:", dim(p), "\n") })
  run("residuals osa", { r <- residuals(f, type = "osa")
    cat("  finite:", all(is.finite(r[, 1])), "\n") })
  run("conditional_effects", { ce <- conditional_effects(f, "x")
    cat("  rows:", nrow(ce[[1]]), "\n") })
  run("emmeans", print(emmeans(f, ~ x)))
  run("summary spec_pars", print(summary(f)$spec_pars))
  run("variables", print(variables(f)))
  run("confint rownames", print(rownames(confint(f))))
  run("bootstrap (default statistic)", {
    bs <- frm_bootstrap(f, nsim = 5, seed = 2)
    cat("  cols:", colnames(bs$t), " NA:", anyNA(bs$t), "\n") })
  run("update() refit", { g <- update(f, data = d[1:250, ])
    cat("  logLik:", as.numeric(logLik(g)), "\n") })
}
# a bootstrap replicate that loses the top category keeps the layout
run("stz: refit on data without the top category", {
  f <- fits$cum_stz
  d2 <- d; d2$y[d2$y == 5L] <- 4L
  g <- frm(y | thres(4) ~ x, data = d2,
           family = cumulative(threshold = "sum_to_zero"))
  cat("  nthres:", family(g)$thres$nthres, "\n")
})
