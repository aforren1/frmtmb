# Reviewer of lane seflags, re-check: the "illcond" reason. Degree-8
# raw polynomial with (1 | g), seeds 3, 8, 10 (the lane's illscan
# fits), and without (1 | g) seed 4: the full warning, the SEs fixef()
# and summary() show, fitted(newdata) and a hypothesis on a coefficient,
# with their warnings, and lmer()/lm() for comparison.
#   Rscript dev/seflags-rev3-ill.R <lane|op2>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev3-common.R")
for (cs in list(c(3, 1), c(8, 1), c(10, 1), c(4, 0))) {
  s <- cs[1]
  re <- cs[2] == 1
  set.seed(s)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
    stats::rnorm(200, 0, 0.2)
  for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
  rhs <- paste0("x", 1:8, collapse = " + ")
  fo <- stats::as.formula(paste("y ~", rhs, if (re) "+ (1 | g)"))
  r <- cap(frm(fo, data = d))
  f <- r$v
  cat(sprintf("\n== degree 8 seed %d re %s | code %d | lost %s\n", s, re,
              f$opt$convergence, lost_of(f)))
  for (w in r$w) cat("  FIT WARNING:", w, "\n")
  fe <- cap(fixef(f))
  cat("  fixef Est.Error:", format(fe$v[, "Est.Error"], digits = 4),
      "| warnings", length(fe$w), "\n")
  sm <- cap(capture.output(print(summary(f))))
  cat(paste0("  | ", head(sm$v, 22)), sep = "\n")
  nd <- d[c(1, 100, 200), ]
  p <- cap(fitted(f, newdata = nd, re_formula = if (re) NA else NULL))
  cat("  fitted(newdata) Est.Error:", format(p$v[, 2], digits = 4),
      "| warnings", length(p$w), "\n")
  for (w in p$w) cat("   PRED WARNING:", w, "\n")
  h <- cap(hypothesis(f, "b_x1 = 0", class = NULL))
  cat("  hypothesis(b_x1 = 0) Est.Error", h$v$hypothesis$Est.Error,
      "| warnings", length(h$w), "\n")
  ref <- if (re) {
    sqrt(diag(as.matrix(stats::vcov(suppressWarnings(suppressMessages(
      lme4::lmer(fo, data = d, REML = FALSE)))))))
  } else sqrt(diag(stats::vcov(stats::lm(fo, data = d))))
  cat("  reference SE:", format(ref, digits = 4), "\n")
}
