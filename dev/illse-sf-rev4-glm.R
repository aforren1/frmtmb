# Lane illse: copy of dev/rel070-seflags-rev4-glm.R by dev/illse-sf-make.R (lane arm: the lane library over
# rellib-r8; base arm: rellib-r8).
# Reviewer of lane seflags, check of punch 2b: the degree-8 raw
# polynomial in non-gaussian families (bernoulli, poisson, cumulative),
# seeds 1 to 10. Per fit: code, lost set, ill warning; per row of
# fitted(newdata) at x = 1.2, 1.5, 1.8, 2, 2.5: the SE against the
# poly(x, 8) fit's (same likelihood, well conditioned; only fits whose
# reference converged without a warning), and whether the call warned.
# A finite row outside 0.9 to 1.1 of the reference with no warning is
# silent and wrong.
#   Rscript dev/seflags-rev4-glm.R <lane|base|op4> <out.tsv>
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-sf-rev4-common.R")
out <- commandArgs(trailingOnly = TRUE)[2]
rhs <- paste0("x", 1:8, collapse = " + ")
nd <- data.frame(x = c(1.2, 1.5, 1.8, 2, 2.5))
for (k in 1:8) nd[[paste0("x", k)]] <- nd$x^k
cat("family\tseed\tcode\tn_lost\till\trow\tse\tref\tcall_warn\tshort\test\test_ref\n", file = out)
for (fam in c("bernoulli", "poisson", "cumulative")) for (s in 1:10) {
  set.seed(s)
  d <- data.frame(x = stats::runif(300, 1, 2))
  eta <- sin(3 * d$x)
  d$y <- switch(fam,
    bernoulli = stats::rbinom(300, 1, stats::plogis(2 * eta)),
    poisson = stats::rpois(300, exp(1 + eta)),
    cumulative = cut(eta + stats::rlogis(300), c(-Inf, -0.5, 0.5, Inf),
                     labels = FALSE))
  for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
  fm <- switch(fam, bernoulli = bernoulli(), poisson = poisson(),
               cumulative = cumulative())
  r <- cap(frm(stats::as.formula(paste("y ~", rhs)), data = d, family = fm))
  g <- cap(frm(y ~ poly(x, 8), data = d, family = fm))
  if (inherits(r$v, "rev_err") || inherits(g$v, "rev_err") ||
      g$v$opt$convergence != 0 || length(g$w)) next
  f <- r$v
  a <- cap(fitted(f, newdata = nd))
  b <- suppressWarnings(fitted(g$v, newdata = nd))
  e <- function(z, col = "Est.Error") if (is.matrix(z)) z[, col] else as.vector(z[, col, ])
  short <- as.numeric(logLik(g$v)) - as.numeric(logLik(f))
  es <- if (inherits(a$v, "rev_err")) rep(NA, length(e(b))) else e(a$v, "Estimate")
  se <- if (inherits(a$v, "rev_err")) rep(NA, length(e(b))) else e(a$v)
  l <- suppressWarnings(ns$sdr_of(f))$se_lost
  for (i in seq_along(se)) {
    cat(fam, s, f$opt$convergence, length(l), any(grepl("ill-conditioned", r$w)),
        i, se[i], e(b)[i], length(a$w), signif(short, 4), es[i],
        e(b, "Estimate")[i], sep = "\t", file = out, append = TRUE)
    cat("\n", file = out, append = TRUE)
  }
}
cat("done\n")
