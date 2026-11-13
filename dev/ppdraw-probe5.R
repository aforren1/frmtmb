# The one hold the core suite made outside this lane's own test:
# test-re-formula-partial.R's fit, y ~ x + (1 + x | g) + (1 | h), seed 3.
# Is theta_2 at a boundary? The slices, the SE check's verdict, and
# lme4's singularity call on the same data.
# Usage: Rscript dev/ppdraw-probe5.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(3)
G <- 15
m <- 8
d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                h = factor(rep(1:5, times = G * m / 5)),
                x = stats::rnorm(G * m))
ug <- stats::rnorm(G, 0, 0.8)
sg <- stats::rnorm(G, 0, 0.4)
uh <- stats::rnorm(5, 0, 0.6)
d$y <- 1 + 0.5 * d$x + ug[d$g] + sg[d$g] * d$x + uh[d$h] +
  stats::rnorm(G * m)
fit <- withCallingHandlers(
  frm(bf(y ~ x + (1 + x | g) + (1 | h)), data = d),
  message = function(m) cat("  M:", conditionMessage(m), "\n"),
  warning = function(w) cat("  W:", conditionMessage(w), "\n"))
p <- fit$opt$par
V <- vcov(fit, full = TRUE)
sd <- sqrt(diag(V))
f0 <- fit$obj$fn(p)
for (j in seq_along(p)) {
  dr <- vapply(c(-1, 1), function(s) {
    q <- p
    q[j] <- p[j] + 2 * s * sd[j]
    fit$obj$fn(q) - f0
  }, 0)
  cat(sprintf("%-28s est %9.4f sd %9.4f drop- %10.4g drop+ %10.4g\n",
              rownames(V)[j], p[j], sd[j], dr[1], dr[2]))
}
cat("se_lost:\n"); print(frmtmb:::sdr_of(fit)$se_lost)
print(VarCorr(fit))
if (requireNamespace("lme4", quietly = TRUE)) {
  lf <- lme4::lmer(y ~ x + (1 + x | g) + (1 | h), data = d, REML = FALSE)
  cat("lme4 isSingular:", lme4::isSingular(lf), "\n")
  print(lme4::VarCorr(lf))
}
