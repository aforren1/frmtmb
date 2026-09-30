# Reviewer re-check: draws_laplace_watch() on hypothesis(), where a later
# draw is non-finite for a reason of its own (data seed 77, draws seed 1).
.libPaths(c("C:/Users/adf44/source/r/wt-sampfix-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
source_lines <- readLines("C:/Users/adf44/source/r/frmtmb-wt-sampfix/dev/sampfix-rev-r2-01-watch.R")
eval(parse(text = source_lines[seq_len(grep("set.seed(77)", source_lines, fixed = TRUE)[1L] - 1L)]))
set.seed(77)
G <- 8; n <- 10
dd <- data.frame(g = factor(rep(seq_len(G), each = n)))
dd$x <- rnorm(nrow(dd))
u <- rnorm(G, 0, 0.6)[dd$g]
dd$y <- 0.5 + 0.4 * dd$x + u + rnorm(nrow(dd), 0, 0.7)
p <- lap_pair(suppressMessages(frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)))
sds <- VarCorr(p$full, summary = FALSE)$g$sd[, 1]
for (h in c(sprintf("exp(1000000000 * (sd_g__Intercept - %.12f)) > 1", sds[1]),
            sprintf("log(%.12f - sd_g__Intercept) > -20", sds[1] + 1e-6))) {
  a <- res(hypothesis(p$lap, h, class = NULL)); b <- res(hypothesis(p$full, h, class = NULL))
  f <- function(r) if (inherits(r, "error")) paste("ERROR:", substr(conditionMessage(r), 1, 100)) else
    sprintf("OK, sample draws non-finite %d of %d", sum(!is.finite(r$samples[[1]])), nrow(r$samples))
  cat(h, "\n  laplace:", f(a), "\n  full:   ", f(b), "\n")
}
