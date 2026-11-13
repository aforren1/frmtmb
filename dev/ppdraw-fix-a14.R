# Fix round: zib zi ~ x (A1) and hurdle_poisson hu ~ x (A4), 4000 draws:
# mean, P(y = 0) and widths against brms 2.23.0's (dev/ppdraw-fix-brms.R:
# A1 5.7036 / 0.0098 / 4.9707 / 7.1353; A4 12.3835 / 0.0033 / 11.6305 /
# 16.7241). Usage: Rscript dev/ppdraw-fix-a14.R <lib or "base">
args <- commandArgs(TRUE)
lib <- if (args[1] != "base") args[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
summ <- function(label, p, fit) {
  w <- function(a, b) mean(apply(p, 2, stats::quantile, b) -
                             apply(p, 2, stats::quantile, a))
  cat(sprintf(paste("%-6s mean %.4f  P(y = 0) %.4f  90%% width %.4f",
                    " 98%% width %.4f  fitted %.4f\n"), label, mean(p),
              mean(p == 0), w(0.05, 0.95), w(0.01, 0.99),
              mean(fitted(fit)[, "Estimate"])))
}
q <- function(expr) suppressMessages(suppressWarnings(expr))
set.seed(7)
dz <- data.frame(x = rnorm(200), N = sample(5:15, 200, TRUE))
dz$y <- rbinom(200, dz$N, plogis(0.3 + 0.5 * dz$x))
f1 <- q(frm(bf(y | trials(N) ~ x, zi ~ x), family = zero_inflated_binomial(),
            data = dz))
set.seed(1)
summ("A1", q(predict(f1, ndraws = 4000, summary = FALSE)), f1)
set.seed(12)
dh <- data.frame(x = rnorm(300))
dh$y <- rpois(300, exp(2.5 + 0.2 * dh$x))
f4 <- q(frm(bf(y ~ x, hu ~ x), family = hurdle_poisson(), data = dh))
set.seed(1)
summ("A4", q(predict(f4, ndraws = 4000, summary = FALSE)), f4)
