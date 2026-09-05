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
L <- diag(rowSums(W)) - W
K <- L + matrix(1 / (1e-3 * n)^2, n, n)
phi <- 1.2 * drop(crossprod(chol(solve(K)), rnorm(n)))
loc <- factor(rep(rownames(W), each = 6), levels = rownames(W))
d <- data.frame(loc = loc, x = rnorm(length(loc)))
d$y <- 1 + 0.5 * d$x + phi[as.integer(d$loc)] + rnorm(nrow(d), 0, 0.5)
fit <- frm(bf(y ~ x + car(W, gr = loc, type = "esicar")) + gaussian(),
           data = d, data2 = list(W = W))

say <- function(lbl, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)),
                warning = function(w) paste("WARNING:", conditionMessage(w)))
  cat(lbl, ": ", if (is.character(r)) r else paste(r, collapse = " "), "\n",
      sep = "")
}
say("importance", tryCatch({
  frm(bf(y ~ x + car(W, gr = loc, type = "esicar")) + gaussian(),
      data = d, data2 = list(W = W), importance = 64); "ACCEPTED"
}, error = function(e) paste("REFUSED:", substr(conditionMessage(e), 1, 90))))

say("predict newdata (all 16 locations)", {
  nd <- data.frame(loc = factor(rownames(W), levels = rownames(W)),
                   x = 0)
  p <- predict(fit, newdata = nd)
  sprintf("n=%d  sum(p - mean(p))=%.3e", length(p), sum(p - mean(p)))
})
say("predict newdata + se", {
  nd <- data.frame(loc = factor(rownames(W)[1:3], levels = rownames(W)), x = 0)
  p <- predict(fit, newdata = nd, se.fit = TRUE)
  paste(round(unlist(lapply(p, function(z) round(as.numeric(z), 6))), 6),
        collapse = " ")
})
say("simulate", {
  set.seed(3); sm <- simulate(fit, nsim = 5)
  sprintf("dim %s  sd %.4f", paste(dim(sm), collapse = "x"),
          stats::sd(as.matrix(sm)))
})
say("frm_sample", {
  set.seed(5)
  sp <- suppressWarnings(frmtmb.sample::frm_sample(
    fit, iter = 400, chains = 1, refresh = 0, seed = 5))
  dr <- as.matrix(posterior::as_draws_matrix(sp))
  bc <- which(startsWith(colnames(dr), "b["))
  fldsum <- apply(dr[, bc, drop = FALSE], 1, sum)
  sprintf("draws %s  b cols %d  max|sum b| %.3e",
          paste(dim(dr), collapse = "x"), length(bc), max(abs(fldsum)))
})
say("ranef condVar", {
  r <- ranef(fit, condVar = TRUE)
  sprintf("sd range %.6f .. %.6f",
          min(attr(r[[1]], "condSD")), max(attr(r[[1]], "condSD")))
})
say("VarCorr", paste(capture.output(print(VarCorr(fit))), collapse = " | "))
say("confint theta", paste(round(as.numeric(confint(fit)["theta_1", ]), 6),
                           collapse = " "))
say("residuals/fitted", sprintf("%.6f %.6f", stats::sd(residuals(fit)),
                                mean(fitted(fit))))

# disconnected graph: per-component hard constraint
W2 <- matrix(0, 12, 12); w1 <- lattice_W(2, 3)
W2[1:6, 1:6] <- w1; W2[7:12, 7:12] <- w1
lv <- paste0("L", 1:12); dimnames(W2) <- list(lv, lv)
set.seed(99)
L2 <- diag(rowSums(W2)) - W2
S <- rbind(c(rep(1, 6), rep(0, 6)), c(rep(0, 6), rep(1, 6)))
K2 <- L2 + t(S) %*% diag(rep(1 / (1e-3 * 6)^2, 2)) %*% S
ph <- drop(crossprod(chol(solve(K2)), rnorm(12)))
lc <- factor(rep(lv, each = 8), levels = lv)
d2 <- data.frame(loc = lc, x = rnorm(length(lc)))
d2$y <- 1 + 0.5 * d2$x + ph[as.integer(d2$loc)] + rnorm(nrow(d2), 0, 0.4)
f2 <- frm(bf(y ~ x + car(W2, gr = loc, type = "esicar")) + gaussian(),
          data = d2, data2 = list(W2 = W2))
re2 <- ranef(f2)[[1]][, 1]
cat(sprintf("disconnected: comp sums %.3e %.3e  logLik %.10f\n",
            sum(re2[1:6]), sum(re2[7:12]), as.numeric(logLik(f2))))
f2i <- frm(bf(y ~ x + car(W2, gr = loc, type = "icar")) + gaussian(),
           data = d2, data2 = list(W2 = W2))
cat(sprintf("disconnected icar logLik %.10f  gap %.4e\n",
            as.numeric(logLik(f2i)),
            as.numeric(logLik(f2)) - as.numeric(logLik(f2i))))
