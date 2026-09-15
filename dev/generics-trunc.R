# Nit 1 of punch round 2. The scale page quoted "a relative 0.1150 at
# the worst row" for the lognormal identity under truncation. Review
# showed that number is max|diff| / max(fitted), that the worst ROW is
# 0.5163, and that the whole thing is an IDENTITY: the truncated mean
# of a lognormal is exp(mu + s^2/2) * Phi(s - a) / Phi(-a) with
# a = (log(lb) - mu) / s, so the per-row relative error of the naive
# formula is 1 - Phi(-a) / Phi(s - a). Reproduced here rather than
# copied, and swept over the bound, because it depends on where the
# bound sits and the page should say that instead of one number.
LIB <- "C:/Users/adf44/source/r/generics-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))

set.seed(2026)
n <- 400
dd <- data.frame(x = rnorm(n), g = factor(rep(1:20, 20)))
eta <- 8 + 0.4 * dd$x + rnorm(20, 0, 0.3)[dd$g]
dd$y <- exp(rnorm(n, eta, 0.4))

one <- function(lb) {
  d3 <- dd[dd$y > lb, ]
  f <- eval(bquote(frm(bf(y | trunc(lb = .(lb)) ~ x + (1 | g)) +
                         lognormal(), data = d3)))
  mu <- predict(f, type = "link")
  s <- sigma(f)
  fv <- fitted(f)
  naive <- exp(mu + s^2 / 2)
  row <- (fv - naive) / fv
  a <- (log(lb) - mu) / s
  closed <- 1 - stats::pnorm(-a) / stats::pnorm(s - a)
  c(lb = lb, rows = nrow(d3),
    old_statistic = max(abs(fv - naive)) / max(abs(fv)),
    worst_row = max(abs(row)), median_row = stats::median(abs(row)),
    a_min = min(a), a_max = max(a),
    identity_resid = max(abs(row - closed)),
    eps_units = max(abs(row - closed)) / .Machine$double.eps)
}
res <- t(vapply(c(500, 2000, 5000), one, numeric(9)))
cat("```\n")
cat("== lognormal under trunc(lb), dev/generics-trunc.R, seed 2026 ==\n")
cat("per-row relative error of exp(mu + s^2/2) against fitted():\n")
cat("  (fitted - naive) / fitted = 1 - Phi(-a) / Phi(s - a),\n")
cat("  a = (log(lb) - mu) / s\n\n")
cat(sprintf("%6s %5s %9s %9s %9s %14s %10s\n", "lb", "rows",
            "old stat", "worst row", "median", "a range",
            "identity"))
for (i in seq_len(nrow(res))) {
  r <- res[i, ]
  cat(sprintf("%6.0f %5.0f %9.4f %9.4f %9.4f %6.2f..%5.2f %10.2e\n",
              r[["lb"]], r[["rows"]], r[["old_statistic"]],
              r[["worst_row"]], r[["median_row"]], r[["a_min"]],
              r[["a_max"]], r[["identity_resid"]]))
}
cat("\n`old stat` is max|diff|/max(fitted), the number the page used to\n")
cat("call the worst row. `identity` is the largest per-row difference\n")
cat("between the measured error and the closed form.\n")
cat("```\n")
