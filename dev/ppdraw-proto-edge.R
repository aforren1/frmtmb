# Prototype of the draw-edge test: along each outer parameter, the
# objective at +-2 marginal standard deviations, others at the estimates.
# The normal draw law predicts a log-likelihood drop of at least 2 there.
# Usage: Rscript dev/ppdraw-proto-edge.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
slice <- function(fit, k = 2) {
  V <- vcov(fit, full = TRUE)
  p <- fit$opt$par
  nm <- frmtmb:::outer_par_names(fit)
  sd <- sqrt(diag(V))[nm]
  f0 <- fit$obj$fn(p)
  out <- t(vapply(seq_along(p), function(j) {
    vapply(c(-1, 1), function(s) {
      q <- p
      q[j] <- p[j] + s * k * sd[j]
      f1 <- tryCatch(fit$obj$fn(q), error = function(e) NA_real_)
      f1 - f0
    }, 0)
  }, c(0, 0)))
  dimnames(out) <- list(nm, c("drop_minus", "drop_plus"))
  cbind(est = p, sd = sd, out)
}
q <- function(expr) suppressMessages(suppressWarnings(expr))

set.seed(7)
n <- 200
dz <- data.frame(x = rnorm(n), N = sample(5:15, n, TRUE))
dz$y <- rbinom(n, dz$N, plogis(0.3 + 0.5 * dz$x))
fz <- q(frm(bf(y | trials(N) ~ x), family = zero_inflated_binomial(),
            data = dz))
cat("== A zib boundary\n"); print(slice(fz))
dz$y2 <- ifelse(runif(n) < 0.3, 0L, dz$y)
fz2 <- q(frm(bf(y2 | trials(N) ~ x), family = zero_inflated_binomial(),
             data = dz))
cat("== A2 zib interior (zi 0.3)\n"); print(slice(fz2))
dz$y3 <- ifelse(runif(n) < 0.02, 0L, dz$y)
fz3 <- q(frm(bf(y3 | trials(N) ~ x), family = zero_inflated_binomial(),
             data = dz))
cat("== A3 zib near boundary (zi 0.02)\n"); print(slice(fz3))

set.seed(7)
do <- data.frame(x = rnorm(400))
ko <- rbinom(400, 1, 0.3)
do$yo <- cut(do$x + ifelse(ko == 1, 2, 0) + rlogis(400),
             c(-Inf, -1, 0.5, 2, Inf), labels = FALSE)
f4 <- q(frm(bf(yo ~ x), family = mixture(cumulative(), cumulative()),
            data = do))
cat("== B ordinal mixture\n"); print(slice(f4))

data(sleepstudy, package = "lme4")
fs <- q(frm(Reaction ~ Days + (Days | Subject), data = sleepstudy))
cat("== C sleepstudy\n")
t0 <- proc.time()[["elapsed"]]
print(slice(fs))
cat("slice time", proc.time()[["elapsed"]] - t0, "\n")
