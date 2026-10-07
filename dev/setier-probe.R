# Lane setier: what the SE check sees on the round's repro fits.
# Usage: Rscript dev/setier-probe.R <lib or "base">
# Run with the reference R and with the OpenBLAS R (dev/setier-openblas.sh)
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
if (identical(args[1], "src")) {
  suppressMessages(pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-setier",
                                     quiet = TRUE))
} else suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), " BLAS probe (s):",
    system.time({m <- matrix(1, 600, 600); m %*% m})[["elapsed"]], "\n")
ns <- asNamespace("frmtmb")
f6 <- function(v) format(signif(v, 6))
conds <- function(expr) {
  w <- character(0); m <- character(0)
  val <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  }, message = function(c) {
    m <<- c(m, conditionMessage(c)); invokeRestart("muffleMessage")
  })
  list(value = val, warnings = w, messages = m)
}
show <- function(tag, fw) {
  fit <- fw$value
  cat("\n==", tag, "code", fit$opt$convergence, "evals", fit$opt$evals,
      "\n")
  for (w in fw$warnings) cat("  W:", substr(w, 1, 300), "\n")
  for (m in fw$messages) cat("  M:", substr(m, 1, 300), "\n")
  cat("  par:", f6(fit$opt$par), "\n")
  h <- fit$cache$hessian_fixed
  if (!is.null(h)) {
    u <- fit$par_units %||% rep(1, nrow(h$H))
    cat("  row max:", f6(apply(abs(h$H * outer(u, u)), 1, max)), "\n")
    cat("  diag:", f6(diag(h$H)), "\n")
    cat("  10*noise:", f6(10 * apply(h$E * outer(u, u), 1, max)), "\n")
  }
  s <- conds(summary(fit))
  for (w in s$warnings) cat("  summary W:", substr(w, 1, 200), "\n")
  sdr <- ns$sdr_of(fit)
  cat("  se:", f6(sqrt(diag(sdr$cov.fixed))), "\n")
  cat("  lost:", paste(names(sdr$se_lost), sdr$se_lost), "\n")
  invisible(fit)
}
data(sleepstudy, package = "lme4")
ss <- sleepstudy
ss$a <- factor(ss$Days %% 3)
show("sleep_nested", conds(frm(Reaction ~ Days + (1 | Subject/a),
                               data = ss)))
show("sleep_slope", conds(frm(Reaction ~ Days + (Days | Subject),
                              data = ss)))
data(cbpp, package = "lme4")
show("cbpp", conds(frm(cbind(incidence, size - incidence) ~ period +
                         (1 | herd), family = binomial(), data = cbpp)))
set.seed(1)
k <- 10
d2 <- data.frame(f = factor(rep(seq_len(k), each = 6)),
                 g = factor(rep(1:6, k)))
d2$y <- rnorm(k)[d2$f] + rnorm(6, 0, 0.5)[d2$g] + rnorm(k * 6, 0, 0.5)
r2 <- show("nl_ridge_re", conds(frm(bf(y ~ a + b, a ~ 0 + f,
                                       b ~ 1 + (1 | g), nl = TRUE),
                                    data = d2)))
p <- conds(frm_linpred(r2, newdata = d2[c(1, 7), ], se.fit = TRUE,
                       dpar = "a"))
cat("  nl_ridge_re pred se:", f6(p$value$se.fit), "warnings",
    length(p$warnings), "\n")
set.seed(514)
d <- data.frame(x = rnorm(240) * 1e-4, z = rnorm(240))
d$yb <- as.integer(d$z > 0)
show("separation", conds(frm(yb ~ z + x, family = bernoulli(), data = d)))
g <- conds(glm(yb ~ z + x, family = binomial(), data = d))
cat("  glm warnings:", g$warnings, " iter", g$value$iter, "\n")
