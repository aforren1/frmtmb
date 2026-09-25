# As phase3b-cleft-probe.R, with the recovery arm's random effects.
# steps, then evaluate the objective and gradient at the start and along
# each coordinate.
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
set.seed(1)
NS <- 8L; NT <- 400L
u <- rnorm(NS, 0, 0.35); b <- rnorm(NS, 0, 0.2)
d <- do.call(rbind, lapply(seq_len(NS), function(s) {
  cond <- rep(0:1, length.out = NT)
  x <- ddm_simulate(NT, mu = 0.4 + 0.9 * cond + u[s],
                    bs = 1.4 * exp(b[s]), ndt = 0.25, bias = 0.5)
  x$cond <- cond; x$s <- factor(s); x
}))
d$code <- 0L; d$y2 <- d$rt
fast <- d$rt < 0.45
d$code[fast] <- -1L; d$rt[fast] <- 0.45
pick <- which(!fast)[seq(5, sum(!fast), by = 5)]
lo <- floor(d$rt[pick] * 10) / 10
d$code[pick] <- 2L; d$y2[pick] <- lo + 0.1; d$rt[pick] <- pmax(lo, 0.45)
fit <- suppressWarnings(frm(bf(rt | dec(upper) + cens(code, y2) ~ cond + (1 | s),
                               bs ~ 1 + (1 | s), ndt ~ 1, bias = 0.5),
                            family = wiener(), data = d, verbose = TRUE,
                            control = frmtmb_control(optCtrl = list(
                              iter.max = 0, eval.max = 1), restarts = 0)))
obj <- fit$obj
p0 <- obj$par
cat("par names:", names(p0), "\n")
cat("start:", format(p0, digits = 4), "\n")
cat("fn:", obj$fn(p0), " gr:", format(obj$gr(p0), digits = 4), "\n")
for (k in seq_along(p0)) for (h in c(-3, -1, -0.3, 0.3, 1, 3, 8)) {
  p <- p0; p[k] <- p[k] + h
  f <- obj$fn(p); g <- obj$gr(p)
  if (!is.finite(f) || any(!is.finite(g))) {
    cat("par", k, names(p0)[k], "step", h, "fn", f, "gr", format(g, digits = 3), "\n")
  }
}
