# Lane osa: trunc() with a random intercept. Under oneStepGeneric the
# truncation normalizer -log Z_j(u) of every later row stays in the
# objective while its density term is switched off. One process of the
# lane build, both routes on the same fit.
#
#   Rscript dev/osa-trunc-re.R
#
# Data seed 31: 12 groups of 8 (the design of test-osa-cdf.R),
# y = 1 + 0.5 x + u_g + e, u_g ~ N(0, 0.8), e ~ N(0, 1), redrawn until
# every row is above 0.5, fitted as gaussian with trunc(lb = 0.5) and
# (1 | g). "exact" is the one-step predictive by quadrature over u (801
# points over +-8 sd) at the fit's estimates: the truncated-normal CDF
# averaged over p(u | earlier rows of the group). The 0.79.1 call
# (oneStepGeneric, range = c(0.5, Inf), splineApprox = FALSE) is made on
# the fit's objective for the first 24 rows only (subset = 1:24, three
# whole groups, each conditioned as in the full run), because each of its
# rows is a quadrature of Laplace evaluations; "cdf" is residuals() on
# every row.
.libPaths(c("C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(31)
G <- 12
K <- 8
n <- G * K
d <- data.frame(g = factor(rep(seq_len(G), each = K)), x = rnorm(n))
u <- rnorm(G, 0, 0.8)[d$g]
d$y <- NA_real_
todo <- seq_len(n)
while (length(todo)) {
  d$y[todo] <- 1 + 0.5 * d$x[todo] + u[todo] + rnorm(length(todo))
  todo <- which(d$y <= 0.5)
}
fit <- frm(y | trunc(lb = 0.5) ~ x + (1 | g), family = gaussian(), data = d)
dp <- lapply(frmtmb:::eval_dpars(fit)[[1]], as.numeric)
b <- fit$estimates$b
eta <- dp$mu - b[as.integer(d$g)]
sig <- dp$sigma[1]
sdu <- exp(fit$estimates$theta[[1]])
ug <- seq(-8, 8, length.out = 801) * sdu
w0 <- dnorm(ug, 0, sdu)
ex <- numeric(n)
for (gi in split(seq_len(n), d$g)) {
  w <- w0
  for (i in gi) {
    m <- eta[i] + ug
    Z <- pnorm((m - 0.5) / sig)
    Fy <- (pnorm((d$y[i] - m) / sig) - pnorm((0.5 - m) / sig)) / Z
    ex[i] <- qnorm(sum(w * Fy) / sum(w))
    w <- w * dnorm(d$y[i], m, sig) / Z
  }
}
cat("BEGIN osa-trunc-re\n")
new <- residuals(fit, type = "osa")[, "Estimate"]
cat(sprintf(paste("cdf          rows 1:%d nan %d max|r - exact| %.3g",
                  "median %.3g mean %.3f sd %.3f\n"),
            n, sum(!is.finite(new)), max(abs(new - ex)),
            median(abs(new - ex)), mean(new), sd(new)))
t0 <- proc.time()[["elapsed"]]
# each row whose integral fails goes through TMB's try(), which prints
# "Error in integrate ... non-finite function value" to try.outFile, and
# becomes NaN; count those lines
tc <- textConnection("trylines", "w")
op <- options(try.outFile = tc)
old <- tryCatch(
  suppressWarnings(RTMB::oneStepPredict(
    fit$obj, observation.name = ".frm_obs", method = "oneStepGeneric",
    discrete = FALSE, range = c(0.5, Inf), splineApprox = FALSE,
    subset = 1:24, trace = FALSE))$residual,
  error = function(e) structure(conditionMessage(e), class = "err"))
options(op)
close(tc)
el <- proc.time()[["elapsed"]] - t0
nfail <- sum(grepl("non-finite function value", trylines))
if (inherits(old, "err")) {
  cat(sprintf(paste("0.79.1 call  rows 1:24 integrate() failed on %d rows;",
                    "oneStepPredict() stopped: %s  %.0fs\n"),
              nfail, old, el))
} else {
  ok <- is.finite(old)
  cat(sprintf(paste("0.79.1 call  rows 1:24 integrate() failed on %d rows,",
                    "nan %d max|r - exact| %.3g median %.3g (finite rows)",
                    "%.0fs\n"),
              nfail, sum(!ok), max(abs(old - ex[1:24])[ok]),
              median(abs(old - ex[1:24])[ok]), el))
}
cat("END osa-trunc-re\n")
