# Review of lane sparsese: fits above options(frmtmb.joint_root_max) that
# are NOT whitened (ML, or a REML design that is not ill-conditioned).
# Prints whether each is whitened, its fitted(re_formula = NA) SEs against
# the GLS conditional value at its variance ratio, and every warning.
#   Rscript sparsese-review-nonwhite.R <groups> <seed>
a <- commandArgs(trailingOnly = TRUE)
ng <- as.integer(a[1])
seed <- as.integer(a[2])
.libPaths(c("C:/Users/adf44/source/r/wt-sparsese-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "groups", ng, "seed", seed, "\n")
set.seed(seed)
n <- 2L * ng
d <- data.frame(x = stats::runif(n, 1, 2),
                g = factor(rep(seq_len(ng), length.out = n)))
d$y <- sin(3 * d$x) + stats::rnorm(ng, 0, 0.5)[d$g] + stats::rnorm(n, 0, 0.3)
for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
nd <- data.frame(x = c(1, 1.5, 2), g = d$g[1])
for (k in 1:7) nd[[paste0("x", k)]] <- nd$x^k
gls_se <- function(fit, deg) {
  of <- fit$opt$par
  s2 <- exp(2 * of[["betad"]])
  tt <- exp(2 * (of[["theta"]] - of[["betad"]]))
  Xo <- if (deg == 0) matrix(1, n, 1) else cbind(1, stats::poly(d$x, deg))
  Xt <- Xo
  for (ix in split(seq_len(n), d$g)) {
    m <- length(ix)
    aa <- (1 - 1 / sqrt(1 + m * tt)) / m
    Xt[ix, ] <- Xo[ix, , drop = FALSE] -
      aa * rep(colSums(Xo[ix, , drop = FALSE]), each = m)
  }
  Ri <- backsolve(qr.R(qr(Xt)), diag(ncol(Xt)))
  Ao <- if (deg == 0) matrix(1, 3, 1) else
    cbind(1, stats::predict(stats::poly(d$x, deg), nd$x))
  sqrt(s2 * rowSums((Ao %*% Ri)^2))
}
cases <- list(
  ml_deg7 = list(fo = y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g),
                 REML = FALSE, deg = 7),
  reml_int = list(fo = y ~ 1 + (1 | g), REML = TRUE, deg = 0),
  reml_lin = list(fo = y ~ x + (1 | g), REML = TRUE, deg = 1))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  w <- character()
  fit <- withCallingHandlers(suppressMessages(frm(cs$fo, data = d,
                                                   REML = cs$REML)),
    warning = function(x) {
      w <<- c(w, substr(conditionMessage(x), 1, 150))
      invokeRestart("muffleWarning")
    })
  sdr <- frmtmb:::sdr_of(fit)
  t0 <- proc.time()[["elapsed"]]
  fp <- withCallingHandlers(fitted(fit, newdata = nd, re_formula = NA),
    warning = function(x) {
      w <<- c(w, substr(conditionMessage(x), 1, 150))
      invokeRestart("muffleWarning")
    })
  el <- proc.time()[["elapsed"]] - t0
  cond <- gls_se(fit, cs$deg)
  cat(sprintf("%-9s whitened %s joint %s lost [%s] fitted %.2f s\n", nm,
              !is.null(fit$beta_whiten),
              NROW(sdr$jointPrecision), paste(sdr$se_lost, collapse = " "),
              el))
  cat("   SE:", format(fp[, "Est.Error"], digits = 8), "\n")
  cat("   SE / GLS conditional - 1:",
      format(fp[, "Est.Error"] / cond - 1, digits = 3), "\n")
  if (length(w)) cat("   W:", paste(unique(w), collapse = " | "), "\n")
}
cat(sprintf("peak working set %.0f MB\n",
            ps::ps_memory_info()[["peak_wset"]] / 2^20))
