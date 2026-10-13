# Reviewer of lane seflags, check of punch 2b: the cumulative fit of
# dev/seflags-rev4-glm.R, seed 9 (on op4 it reaches the maximum): the
# category probabilities' SEs (fitted(), the finite-difference route)
# against the linear predictor's SE (frm_linpred(), the analytic route),
# both against the poly(x, 8) fit.
#   Rscript dev/seflags-rev4-cum.R <lane|op4|base>
source("C:/Users/adf44/source/r/frmtmb-wt-release/dev/rel070-seflags-rev4-common.R")
rhs <- paste0("x", 1:8, collapse = " + ")
nd <- data.frame(x = c(1.2, 1.5, 1.8, 2, 2.5))
for (k in 1:8) nd[[paste0("x", k)]] <- nd$x^k
for (s in c(9, 2, 5)) {
  set.seed(s)
  d <- data.frame(x = stats::runif(300, 1, 2))
  eta <- sin(3 * d$x)
  d$y <- cut(eta + stats::rlogis(300), c(-Inf, -0.5, 0.5, Inf), labels = FALSE)
  for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
  r <- cap(frm(stats::as.formula(paste("y ~", rhs)), data = d, family = cumulative()))
  f <- r$v
  g <- suppressWarnings(frm(y ~ poly(x, 8), data = d, family = cumulative()))
  cat(sprintf("\n== seed %d | code %d | lost %s | short %.4g\n", s, f$opt$convergence,
              lost_of(f), as.numeric(logLik(g)) - as.numeric(logLik(f))))
  for (w in r$w) cat("  W:", substr(w, 1, 160), "\n")
  a <- cap(frm_linpred(f, newdata = nd, se.fit = TRUE))
  b <- suppressWarnings(frm_linpred(g, newdata = nd, se.fit = TRUE))
  cat("  eta SE / ref     :", format(a$v$se.fit / b$se.fit, digits = 4), "| warn", length(a$w), "\n")
  a <- cap(fitted(f, newdata = nd))
  b <- suppressWarnings(fitted(g, newdata = nd))
  cat("  P(y = 1) SE / ref:", format(a$v[, "Est.Error", 1] / b[, "Est.Error", 1], digits = 4),
      "| warn", length(a$w), "\n")
  cat("  P(y = 1) est / ref:", format(a$v[, "Estimate", 1] / b[, "Estimate", 1], digits = 4), "\n")
}
