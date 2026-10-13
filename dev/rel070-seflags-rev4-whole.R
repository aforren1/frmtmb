# Reviewer of lane seflags, check of punch 2b: the consumers that read
# the covariance WHOLE (vcov(), get_vcov() for marginaleffects, the
# emmeans bridge), on ill-conditioned fits without random effects
# (degree-8 raw polynomial through I(x^k), seeds 1 to 12, the code-0
# fits). Against the same call on the poly(x, 8) fit. Per fit: the ill
# warning, and per route the ratio at x = 1, 1.5, 2, 2.5 and whether
# frmtmb warned.
#   Rscript dev/seflags-rev4-whole.R <lane|base|op4>
source("C:/Users/adf44/source/r/frmtmb-wt-release/dev/rel070-seflags-rev4-common.R")
suppressMessages({library(emmeans); library(marginaleffects)})
pw <- paste0(" + I(x^", 2:8, ")", collapse = "")
nd <- data.frame(x = c(1, 1.5, 2, 2.5))
fw <- function(w) sum(!grepl("not known to be supported", w))
for (s in 1:12) {
  set.seed(s)
  d <- data.frame(x = stats::runif(200, 1, 2))
  d$y <- sin(3 * d$x) + stats::rnorm(200, 0, 0.3)
  r <- cap(frm(stats::as.formula(paste0("y ~ x", pw)), data = d))
  f <- r$v
  g <- suppressWarnings(frm(y ~ poly(x, 8), data = d))
  ill <- any(grepl("ill-conditioned", r$w))
  cat(sprintf("\n== seed %2d | code %d | lost %s | ill warning %s | other fit warnings %d\n",
              s, f$opt$convergence, lost_of(f), ill, length(r$w) - ill))
  if (f$opt$convergence != 0) next
  q <- function(expr, lab) {
    a <- cap(eval(expr, list(m = f)))
    b <- cap(eval(expr, list(m = g)))
    v <- if (inherits(a$v, "rev_err") || inherits(b$v, "rev_err")) {
      paste("ERR", substr(c(a$v, b$v)[1], 1, 80))
    } else paste(format(a$v / b$v, digits = 4), collapse = " ")
    cat(sprintf("  %-30s %s | frmtmb warnings %d\n", lab, v, fw(a$w)))
  }
  q(quote(fitted(m, newdata = nd)[, "Est.Error"]), "fitted()")
  q(quote({
    b <- frm_lp_basis(m, newdata = nd)
    A <- as.matrix(b$A)
    V <- vcov(m)
    sqrt(pmax(rowSums((A %*% V) * A), 0))
  }), "A vcov() A'")
  q(quote(marginaleffects::predictions(m, newdata = nd)$std.error),
    "marginaleffects predictions()")
  q(quote(marginaleffects::slopes(m, variables = "x", newdata = nd)$std.error),
    "marginaleffects slopes(x)")
  q(quote(marginaleffects::comparisons(m, variables = list(x = 0.5),
                                       newdata = nd)$std.error),
    "marginaleffects comparisons(+0.5)")
  q(quote(summary(emmeans::emmeans(m, ~ x, at = list(x = nd$x)))$SE),
    "emmeans()")
  q(quote(summary(emmeans::emtrends(m, ~ x, var = "x",
                                    at = list(x = nd$x)))$SE),
    "emtrends(x)")
}
