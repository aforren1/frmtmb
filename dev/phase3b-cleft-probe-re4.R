# As phase3b-cleft-probe-re.R, on the unfitted objective (dry_run).
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
                            family = wiener(), data = d, dry_run = "objective",
                            control = frmtmb_control(optCtrl = list(
                              iter.max = 0, eval.max = 1), restarts = 0)))
fam <- frmtmb::single_response(fit, "fit")$family
av <- fit$frame$aterm_values[["rt"]]
y <- as.numeric(fit$frame$y[["rt"]])
for (mu in c(0, 5, 20, -20)) {
  dp <- list(mu = mu, bs = 1.5, ndt = 0.225, bias = 0.5)
  lf <- fam$lpdf(y, dp, av)
  Fv <- fam$lcdf(y, dp, av)
  F2 <- fam$lcdf(av$cens_y2, dp, av)
  cen <- av$cens
  cat("mu", mu, ": dens non-finite", sum(!is.finite(lf[cen == 0])),
      " left log(Fv) non-finite", sum(!is.finite(log(Fv[cen == -1]))),
      " interval log(F2-Fv) non-finite", sum(!is.finite(log(F2[cen == 2] - Fv[cen == 2]))),
      " any NaN in Fv", anyNA(Fv), " F2", anyNA(F2), "\n")
}
cat("names(av):", names(av), "\n")
cat("identical y2 test:", identical(as.numeric(av$cens_y2), as.numeric(av$cens_y2)), "\n")
