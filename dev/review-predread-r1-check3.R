# Short check, item 3: marginaleffects on both builds and on brms, seed 7402 (me), 7403 (gp).
# collapse, which brms's route needs, is in rev-predread-lib/extra, the review's own.
#   REV_LIB=before Rscript dev/review-predread-r1-check3.R
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
.libPaths(c(.libPaths()[1L], "C:/Users/adf44/source/r/rev-predread-lib/extra", .libPaths()[-1L]))
set.seed(7402)
n <- 200
d <- data.frame(z = rnorm(n), lat = rnorm(n), sx = 0.4, f = factor(sample(c("a", "b"), n, TRUE)))
d$xo <- d$lat + rnorm(n, 0, 0.4)
d$y <- 1 + 1.2 * d$lat - 0.5 * d$z + 0.6 * (d$f == "b") + rnorm(n, 0, 0.5)
P <- rev_pair("bf(y ~ me(xo, sx) + z + f)", "gaussian()", d, ndraws = 4000)
fit <- P$fit
fx <- P$fx
print(signif(summary(fit)$fixed[, 1:2], 4))
ndm <- data.frame(xo = c(-1, 1), sx = 0.6, z = 0, f = "a")
show <- function(lab, e) {
  r <- tryCatch(suppressWarnings(suppressMessages(e)), error = function(err) substr(conditionMessage(err), 1, 100))
  if (is.character(r)) {
    cat(sprintf("%-26s %s\n", lab, r))
  } else {
    cat(sprintf("%-26s est %s  se %s  ci %s %s\n", lab, paste(signif(r$estimate, 4), collapse = " "),
                paste(signif(as.numeric(r$std.error %||% NA), 4), collapse = " "), paste(signif(r$conf.low, 4), collapse = " "),
                paste(signif(r$conf.high, 4), collapse = " ")))
  }
}
cat("--- frmtmb\n")
show("predictions(newdata)", marginaleffects::predictions(fit, newdata = ndm))
show("avg_predictions", marginaleffects::avg_predictions(fit))
show("avg_slopes(z)", marginaleffects::avg_slopes(fit, variables = "z"))
show("avg_slopes(xo)", marginaleffects::avg_slopes(fit, variables = "xo"))
show("avg_comparisons(f)", marginaleffects::avg_comparisons(fit, variables = "f"))
show("hypotheses", marginaleffects::hypotheses(fit, hypothesis = "b2 = 0"))
cat("--- brms at fixed values, no coefficient uncertainty\n")
if (!identical(Sys.getenv("REV_LIB"), "before")) {
  show("predictions(newdata)", marginaleffects::predictions(fx, newdata = ndm))
  show("avg_predictions", marginaleffects::avg_predictions(fx))
  show("avg_slopes(z)", marginaleffects::avg_slopes(fx, variables = "z"))
  show("avg_slopes(xo)", marginaleffects::avg_slopes(fx, variables = "xo"))
  show("avg_comparisons(f)", marginaleffects::avg_comparisons(fx, variables = "f"))
}
cat("\n--- exact gp(): what the delta method can move\n")
set.seed(7403)
m <- 80
g <- data.frame(x = sort(runif(m, 0, 6)), z = rnorm(80))
g$y <- sin(g$x) + 0.3 * g$z + rnorm(m, 0, 0.25)
nd <- data.frame(x = c(g$x[5], 3.33, 7, 9), z = 0)
for (sp in list(bf(y ~ z + gp(x)), bf(y ~ z + gp(x, k = 15, c = 5/4)), bf(y ~ z + s(x, k = 8)), bf(y ~ z + x + (1 | z2)))) {
  dd <- g
  dd$z2 <- factor(rep(1:8, 10))
  nd$z2 <- factor(1, levels = 1:8)
  f1 <- suppressWarnings(suppressMessages(frm(sp, family = gaussian(), data = dd)))
  cat(deparse(sp$formula), "\n")
  cat("  get_coef length", length(marginaleffects::get_coef(f1)), " names", paste(utils::head(names(marginaleffects::get_coef(f1)), 5), collapse = " "), "\n")
  l <- frm_linpred(f1, newdata = nd, se.fit = TRUE)
  cat("  frm_linpred se", signif(l$se.fit, 4), "\n")
  r <- tryCatch(suppressWarnings(suppressMessages(marginaleffects::predictions(f1, newdata = nd))), error = function(e) conditionMessage(e))
  cat("  predictions se", if (is.character(r)) substr(r, 1, 90) else signif(r$std.error, 4), "\n")
}
