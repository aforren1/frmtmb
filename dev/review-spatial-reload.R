# Reviewer: a saved car() fit read back without frmtmb.spatial.
# Usage: Rscript review-spatial-reload.R <mode fit|reload> <libs>
a <- commandArgs(TRUE)
.libPaths(c(strsplit(a[2], ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
out <- "C:/Users/adf44/source/r/frmtmb-wt-spatial/dev/review-spatial-log/"
suppressMessages(library(frmtmb))
cat("lib:", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), " spatial installed:",
    nzchar(system.file(package = "frmtmb.spatial")), "\n")
lattice_W <- function(r, c) {
  g <- expand.grid(r = seq_len(r), c = seq_len(c))
  n <- nrow(g)
  W <- 1 * (abs(outer(g$r, g$r, "-")) + abs(outer(g$c, g$c, "-")) == 1)
  dimnames(W) <- list(paste0("L", seq_len(n)), paste0("L", seq_len(n)))
  W
}
if (a[1] == "fit") {
  suppressMessages(library(frmtmb.spatial))
  set.seed(42)
  W <- lattice_W(4, 4)
  loc <- factor(rep(rownames(W), each = 6), levels = rownames(W))
  d <- data.frame(loc = loc, x = rnorm(length(loc)))
  d$y <- 1 + 0.5 * d$x + rnorm(16)[as.integer(loc)] + rnorm(nrow(d), 0, 0.5)
  fit <- frm(bf(y ~ x + car(W, gr = loc, type = "escar")) + gaussian(),
             data = d, data2 = list(W = W))
  saveRDS(fit, paste0(out, "car-fit.rds"))
  cat("saved; logLik", format(as.numeric(logLik(fit)), digits = 10), "\n")
  # one row per area, gaussian: the OLRE check and a coupled field
  W6 <- lattice_W(6, 6)
  d1 <- data.frame(loc = factor(rownames(W6), levels = rownames(W6)),
                   x = rnorm(36))
  Q <- diag(rowSums(W6)) - 0.9 * W6
  phi <- drop(crossprod(chol(solve(Q)), rnorm(36)))
  d1$y <- 0.5 * d1$x + phi + rnorm(36, 0, 0.5)
  for (ty in c("escar", "icar", "bym2")) {
    r <- tryCatch(withCallingHandlers({
      f1 <- frm(bf(y ~ x + car(W6, gr = loc, type = ty)) + gaussian(),
                data = d1, data2 = list(W6 = W6))
      paste("FIT logLik", format(as.numeric(logLik(f1)), digits = 10))
    }, warning = function(w) {
      cat("  [warning]", ty, substr(conditionMessage(w), 1, 220), "\n")
      invokeRestart("muffleWarning")
    }, message = function(m) {
      cat("  [message]", ty, substr(conditionMessage(m), 1, 220), "\n")
      invokeRestart("muffleMessage")
    }), error = function(e) paste("ERR", conditionMessage(e)))
    cat("one row per area, type", ty, "->", r, "\n")
  }
  # gr(cov = A) control on the same one-row-per-level design
  A <- solve(Q); dimnames(A) <- dimnames(W6)
  r <- tryCatch(withCallingHandlers({
    f2 <- frm(y ~ x + (1 | gr(loc, cov = A)), data = d1,
              data2 = list(A = A))
    paste("FIT logLik", format(as.numeric(logLik(f2)), digits = 10))
  }, warning = function(w) {
    cat("  [warning] gr_cov", substr(conditionMessage(w), 1, 220), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) paste("ERR", conditionMessage(e)))
  cat("one row per level, gr(cov = A) ->", r, "\n")
} else {
  fit <- readRDS(paste0(out, "car-fit.rds"))
  cat("after readRDS, frmtmb.spatial loaded:",
      isNamespaceLoaded("frmtmb.spatial"), "\n")
  for (m in c("print", "summary", "logLik", "VarCorr", "vcov", "fitted",
              "predict_new", "simulate", "ranef", "confint_varcorr",
              "update")) {
    r <- tryCatch({
      switch(m,
             print = capture.output(print(fit)),
             summary = capture.output(summary(fit)),
             logLik = logLik(fit),
             VarCorr = VarCorr(fit),
             vcov = vcov(fit),
             fitted = fitted(fit),
             predict_new = predict(fit, newdata = data.frame(
               x = 0, loc = "Lnew"), allow_new_levels = TRUE),
             simulate = simulate(fit, nsim = 1, seed = 1),
             ranef = ranef(fit),
             confint_varcorr = confint_varcorr(fit),
             update = update(fit, . ~ . - x))
      "ok"
    }, error = function(e) paste("ERR", substr(conditionMessage(e), 1, 250)))
    cat(sprintf("%-16s %s\n", m, r))
  }
}
