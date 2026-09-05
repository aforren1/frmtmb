# P2: con_sd does not touch the esicar LIKELIHOOD but does reach every
# delta-method standard error, because lp_delta_A() uses dc/db = I.
# Measure the size so the doc sentences and the test can state it.
lib_es <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/es-lib"
.libPaths(c(lib_es, "C:/Users/adf44/AppData/Local/R/win-library/4.6", .libPaths()))
library(frmtmb)
lattice_W <- function(r, c) {
  g <- expand.grid(r = seq_len(r), c = seq_len(c)); n <- nrow(g)
  W <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1) W[i, j] <- 1
  }
  dimnames(W) <- list(paste0("L", seq_len(n)), paste0("L", seq_len(n)))
  W
}
set.seed(42); W <- lattice_W(4, 4); n <- nrow(W)
K <- diag(rowSums(W)) - W + matrix(1 / (1e-3 * n)^2, n, n)
phi <- 1.2 * drop(crossprod(chol(solve(K)), rnorm(n)))
loc <- factor(rep(rownames(W), each = 6), levels = rownames(W))
d <- data.frame(loc = loc, x = rnorm(length(loc)))
d$y <- 1 + 0.5 * d$x + phi[as.integer(d$loc)] + rnorm(nrow(d), 0, 0.5)
f <- function(cs) frm(bf(y ~ x + car(W, gr = loc, type = "esicar",
                                     con_sd = cs)) + gaussian(),
                      data = d, data2 = list(W = W))
se <- function(fit) predict(fit, se.fit = TRUE)$se.fit
cs <- c(1e-2, 1e-3, 1e-4)
fits <- lapply(cs, f)
ses <- lapply(fits, se)
cat("logLik (must not move):\n")
for (i in seq_along(cs)) {
  cat(sprintf("  %.0e  %.12f\n", cs[i], as.numeric(logLik(fits[[i]]))))
}
cat("\nse.fit range and the excess over the 1e-4 fit:\n")
for (i in seq_along(cs)) {
  cat(sprintf("  %.0e  %.9f .. %.9f\n", cs[i], min(ses[[i]]), max(ses[[i]])))
}
cat("\nvariance excess, se(cs)^2 - se(1e-4)^2, against cs^2 - 1e-8:\n")
for (i in seq_along(cs)) {
  ex <- ses[[i]]^2 - ses[[3]]^2
  cat(sprintf("  %.0e  measured %.6e  predicted %.6e  max dev %.2e\n",
              cs[i], mean(ex), cs[i]^2 - 1e-8, max(abs(ex - (cs[i]^2 - 1e-8)))))
}
cat("\nrelative SE inflation over the exact Jacobian, se/sqrt(se^2 - cs^2) - 1:\n")
for (i in seq_along(cs)) {
  rel <- ses[[i]] / sqrt(ses[[i]]^2 - cs[i]^2) - 1
  cat(sprintf("  %.0e  %.6e .. %.6e\n", cs[i], min(rel), max(rel)))
}
cat("\nranef condSD range:\n")
for (i in seq_along(cs)) {
  sd_i <- attr(ranef(fits[[i]], condVar = TRUE)[[1]], "condSD")
  cat(sprintf("  %.0e  %.9f .. %.9f\n", cs[i], min(sd_i), max(sd_i)))
}
cat("\nVarCorr sd(car) vs exp(theta1) (must be clean):\n")
for (i in seq_along(cs)) {
  v <- sqrt(unname(VarCorr(fits[[i]])[[1]][1, 1]))
  cat(sprintf("  %.0e  %.15f  %.3e\n", cs[i], v,
              v - exp(fits[[i]]$estimates$theta[[1]])))
}
cat("\nwhat con_sd = 0.1 would do:\n")
f01 <- f(0.1)
s01 <- se(f01)
cat(sprintf("  se range %.6f .. %.6f, logLik %.12f\n", min(s01), max(s01),
            as.numeric(logLik(f01))))
cat(sprintf("  inflation over the 1e-4 fit: %.4f%% .. %.4f%%\n",
            100 * (min(s01 / ses[[3]]) - 1), 100 * (max(s01 / ses[[3]]) - 1)))
