# Release review 0.80.0: readers' last change (irrational steps, the net
# row by row). One row per fitted() call, so a warning is that row's own.
# sin(a)^2 + c at a = 0 with a draw sd s on a: exact mean
# (1 - exp(-2 s^2)) / 2 + c. Also a lognormal-free smooth row and the
# 10 pi row.
#
#   Rscript dev/review-rel-readers.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(strsplit(lib, ";")[[1]], "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "\n")
set.seed(7208)
n <- 800
d <- data.frame(l1 = rnorm(n), l2 = rnorm(n), s1 = 0.2, s2 = 0.2)
d$x1 <- d$l1 + rnorm(n, 0, 0.2)
d$x2 <- d$l2 + rnorm(n, 0, 0.2)
d$y <- sin(0.3 + 0.8 * d$l1)^2 + (0.5 + 0.5 * d$l2) + rnorm(n, 0, 0.1)
fit <- frm(bf(y ~ sin(a)^2 + c, a ~ me(x1, s1), c ~ me(x2, s2), nl = TRUE),
           family = gaussian(), data = d)
lp <- function(nd, dp) as.numeric(frm_linpred(fit, newdata = nd, dpar = dp))
q <- data.frame(x1 = 0, x2 = 0, s1 = 1, s2 = 1)
ba <- lp(transform(q, x1 = 1), "a") - lp(q, "a")
x0 <- -lp(q, "a") / ba
for (s in c(0.3, 1, 2, pi / 2, pi, 1.5 * pi, 2 * pi, 4, 10 * pi)) {
  nd <- data.frame(x1 = x0, x2 = 0, s1 = s / abs(ba), s2 = 0.5)
  ws <- character()
  f <- withCallingHandlers(fitted(fit, newdata = nd), warning = function(w) {
    ws <<- c(ws, class(w)[1L])
    invokeRestart("muffleWarning")
  })
  exact <- (1 - exp(-2 * s^2)) / 2 + lp(nd, "c")
  cat(sprintf("sd %7.4f fitted %.6f exact %.6f rel %+.2e warn %s\n", s,
              f[, "Estimate"], exact, f[, "Estimate"] / exact - 1,
              if (length(ws)) paste(ws, collapse = ",") else "-"))
}
# two rows in one call: the settled row's value must not depend on the
# unsettled row beside it
nd2 <- data.frame(x1 = x0, x2 = 0, s1 = c(pi, 10 * pi) / abs(ba), s2 = 0.5)
f2 <- suppressWarnings(fitted(fit, newdata = nd2))
nd1 <- nd2[1, ]
f1 <- fitted(fit, newdata = nd1)
cat(sprintf("pi row alone %.8f / %.6f, beside 10 pi %.8f / %.6f\n",
            f1[1, "Estimate"], f1[1, "Est.Error"], f2[1, "Estimate"],
            f2[1, "Est.Error"]))
cat("READERS DONE\n")
