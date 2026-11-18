# Lane sparsese: item 4952, after dev/reviews/selimit2-short-3370ref.R.
# A gaussian REML fit, degree-7 raw polynomial, (1 | g) with `ng` groups of
# 2 rows. Prints the conditional SE of fitted(re_formula = NA) at x = 1,
# 1.5, 2 (GLS at the fitted group variance in an orthogonal basis; a fitted
# SE adds the variance parameters' propagation and must not fall below it),
# the fitted SEs, and the time and memory of the first fitted() call and of
# vcov().
#   Rscript sparsese-3370ref.R <base|lane> <groups> <seed> [limit]
# `limit`, when given, is set as options(frmtmb.joint_root_max) before the
# fit, so that it decides the route of the fit's factor and of fitted().
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
ng <- as.integer(a[2])
seed <- as.integer(a[3])
lim <- if (length(a) >= 4L) as.numeric(a[4]) else NA
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-sparsese-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "lib", dirname(find.package("frmtmb")), "groups", ng,
    "seed", seed, "limit", lim, "\n")
peak <- function() ps::ps_memory_info()[["peak_wset"]] / 2^20
wset <- function() ps::ps_memory_info()[["wset"]] / 2^20
set.seed(seed)
n <- 2L * ng
d <- data.frame(x = stats::runif(n, 1, 2),
                g = factor(rep(seq_len(ng), length.out = n)),
                f = factor(rep(c("a", "b"), length.out = n)))
d$y <- sin(3 * d$x) + stats::rnorm(n, 0, 0.3)
for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
if (!is.na(lim)) options(frmtmb.joint_root_max = lim)
t0 <- proc.time()[["elapsed"]]
fit <- withCallingHandlers(suppressMessages(frm(
  y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g), data = d, REML = TRUE)),
  warning = function(w) {
    cat("fit W:", substr(conditionMessage(w), 1, 200), "\n")
    invokeRestart("muffleWarning")
  })
cat(sprintf("fit: %.2f s, peak working set %.0f MB\n",
            proc.time()[["elapsed"]] - t0, peak()))
cat("joint parameters:", nrow(frmtmb:::sdr_of(fit)$jointPrecision), "\n")
est <- fit$opt$par
cat("par:", paste(names(est), format(est, digits = 6)), "\n")
lost <- frmtmb:::sdr_of(fit)$se_lost
cat("lost:", paste(names(lost), lost), "\n")
cat("joint_root:", if (is.null(frmtmb:::sdr_of(fit)$joint_root)) "NULL" else
  paste(dim(frmtmb:::sdr_of(fit)$joint_root), collapse = " x "), "\n")
s2 <- exp(2 * fit$estimates$betad[1])
t <- exp(2 * (fit$estimates$theta[1] - fit$estimates$betad[1]))
pb <- stats::poly(d$x, 7)
X <- cbind(1, pb)
nd <- data.frame(x = c(1, 1.5, 2))
A <- cbind(1, stats::predict(pb, nd$x))
# whitened GLS: subtract a per group, as test-se-reparam-reml.R's reml_gls
Xt <- X
for (ix in split(seq_len(n), d$g)) {
  m <- length(ix)
  aa <- (1 - 1 / sqrt(1 + m * t)) / m
  Xt[ix, ] <- X[ix, , drop = FALSE] - aa * rep(colSums(X[ix, , drop = FALSE]),
                                               each = m)
}
Ri <- backsolve(qr.R(qr(Xt)), diag(ncol(Xt)))
cond <- sqrt(s2 * rowSums((A %*% Ri)^2))
cat("variance ratio t", format(t, digits = 4), "\n")
cat("conditional SE:", format(cond, digits = 10), "\n")
nd2 <- data.frame(x = c(1, 1.5, 2), g = d$g[1])
for (k in 1:7) nd2[[paste0("x", k)]] <- nd2$x^k
invisible(gc(reset = TRUE))
w0 <- wset()
p0 <- peak()
t0 <- proc.time()[["elapsed"]]
fp <- withCallingHandlers(stats::fitted(fit, newdata = nd2, re_formula = NA),
  warning = function(w) {
    cat("fitted W:", substr(conditionMessage(w), 1, 200), "\n")
    invokeRestart("muffleWarning")
  })
el <- proc.time()[["elapsed"]] - t0
g <- gc()
cat(sprintf(paste0("fitted(): %.2f s; R heap max used %.0f MB; working set ",
                   "before %.0f MB, peak before %.0f MB, peak after %.0f MB\n"),
            el, sum(g[, ncol(g)]), w0, p0, peak()))
se <- fp[, "Est.Error"]
cat("fitted SE:", format(se, digits = 10), "\n")
cat("SE / conditional - 1:", format(se / cond - 1, digits = 4), "\n")
t0 <- proc.time()[["elapsed"]]
V <- withCallingHandlers(stats::vcov(fit),
  warning = function(w) {
    cat("vcov W:", substr(conditionMessage(w), 1, 200), "\n")
    invokeRestart("muffleWarning")
  })
cat(sprintf("vcov(): %.2f s, peak %.0f MB\n",
            proc.time()[["elapsed"]] - t0, peak()))
# a' V a from vcov() in the raw basis, against the same rows
Araw <- cbind(1, sapply(1:7, function(k) nd$x^k))
cat("vcov SE (a' V a):",
    format(sqrt(pmax(rowSums((Araw %*% V) * Araw), 0)), digits = 10), "\n")
cat("vcov SE:", format(sqrt(diag(V)), digits = 10), "\n")
