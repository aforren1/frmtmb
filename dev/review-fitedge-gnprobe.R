# Review of lane fitedge, item 1: the Gauss-Newton linearization where the
# start is far off, and where the body is flat at the start in a coefficient.
# Args: arm (base | lane). Data: dev/review-fitedge-designs.R, seed 1, n 60.
# For each case: the default fit, and the best of 12 theta starts (the same
# user start of beta, so the reference answers the same question).
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
source(file.path(FE_WT, "dev/review-fitedge-designs.R"))
ns <- asNamespace("frmtmb")
one <- function(lab, fml, d, fam, st = NULL, nth = 2L) {
  t0 <- proc.time()[[3]]
  r <- if (is.null(st)) fe_fit(fml, data = d, family = fam) else
    fe_fit(fml, data = d, family = fam, start = st)
  tt <- proc.time()[[3]] - t0
  ll <- if (inherits(r$fit, "error")) NA else as.numeric(logLik(r$fit))
  best <- -Inf
  for (rho in c(0.15, 0.5, 1.5, 5)) for (sd0 in c(0.3, 1, 3)) {
    s2 <- c(st, list(theta = rep(c(log(sd0), log(rho)), length.out = nth)))
    q <- fe_fit(fml, data = d, family = fam, start = s2)
    if (!inherits(q$fit, "error")) best <- max(best, as.numeric(logLik(q$fit)))
  }
  cat(sprintf("GN %s %s | ll %.4f ref %.4f short %.4f | conv %s nwarn %d | %.2f s%s\n",
              arm, lab, ll, best, ll - best,
              if (inherits(r$fit, "error")) "ERR" else r$fit$opt$convergence,
              length(r$warn), tt,
              if (inherits(r$fit, "error")) paste(" ERROR",
                conditionMessage(r$fit)) else ""))
  for (w in r$warn) cat("    W:", substr(w, 1, 120), "\n")
  # what the body start gives on the lane, and from which template
  if (arm == "lane" && !inherits(r$fit, "error")) {
    fr <- r$fit$frame
    tpl <- fr$par_template
    if (!is.null(st$beta)) tpl$beta[names(st$beta)] <- st$beta
    w <- tryCatch(ns$gp_start_working(fr, tpl), error = function(e) e)
    cat("    body start:", if (inherits(w, "error")) conditionMessage(w) else
      if (is.null(w)) "NULL" else paste(round(w$theta, 3), collapse = " "),
      "| fitted theta", paste(round(r$fit$opt$par[names(r$fit$opt$par) ==
                                                     "theta"], 3),
                              collapse = " "), "\n")
  }
}
d <- rv_mk("E", 1, 60, 2)
fE <- rv_spell$E$nl
one("E default", fE, d, gaussian())
# far off: b = 3 makes exp(b z) reach e^9 at z = 3
one("E start b=3", fE, d, gaussian(), list(beta = c("b_(Intercept)" = 3)))
one("E start b=-3", fE, d, gaussian(), list(beta = c("b_(Intercept)" = -3)))
one("E start a=50", fE, d, gaussian(), list(beta = c("a_(Intercept)" = 50)))
# flat: a * exp(b z) with b ~ gp(x) and a at 0, so the body does not read b
dL <- rv_mk("L", 1, 60, 2)
one("L default (a starts 0)", rv_spell$L$nl, dL, gaussian())
one("L start a=2", rv_spell$L$nl, dL, gaussian(),
    list(beta = c("a_(Intercept)" = 2)))
# flat in the field's own parameter: a * b with a ~ gp(x), b at 0
dM <- rv_mk("G", 1, 60, 2)
fM <- bf(y ~ b * a + c * z^2, a ~ gp(x), b ~ 1, c ~ 1, nl = TRUE)
one("M b*a default (b starts 0)", fM, dM, gaussian())
one("M b*a start b=1", fM, dM, gaussian(), list(beta = c("b_(Intercept)" = 1)))
# a field where the body saturates: a / (1 + exp(-b (z - c))) with c far off
dF <- rv_mk("F", 1, 60, 2)
one("F start c=8", rv_spell$F$nl, dF, gaussian(),
    list(beta = c("c_(Intercept)" = 8)))
one("F start b=20", rv_spell$F$nl, dF, gaussian(),
    list(beta = c("b_(Intercept)" = 20)))
# y on a large scale: exp(a) with a near 9
dC <- rv_mk("E", 1, 60, 1)
dC$y <- 5000 * dC$y
one("E y*5000", fE, dC, gaussian())
cat("DONE\n")
