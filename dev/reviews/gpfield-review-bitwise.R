# Review of lane gpfield: exp_quad fits lane against 0.82.0, bit for bit
# (estimates, logLik, SEs, warnings), each arm in its own process.
# Usage: Rscript gpfield-review-bitwise.R <lane|base> <tag>
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
a <- commandArgs(TRUE)
arm <- a[1]; tag <- a[2]
gf_libs(arm)
mk <- function(seed, n) {
  set.seed(seed)
  x <- sort(round(runif(n, 0, 10), 3)); x2 <- round(runif(n, 0, 5), 2)
  g <- factor(sample(c("a", "b"), n, TRUE))
  f <- sin(x)
  data.frame(x, x2, g, yg = f + 0.3 * x2 + 0.25 * rnorm(n),
             yp = rpois(n, exp(1.2 + 0.6 * f)),
             yb = rbinom(n, 1, plogis(1.5 * f)))
}
d <- mk(9500, 60)
rows <- list(
  g1 = list(bf(yg ~ gp(x)), gaussian()),
  g2ani = list(bf(yg ~ gp(x, x2, iso = FALSE)), gaussian()),
  p1 = list(bf(yp ~ gp(x)), poisson()),
  b1 = list(bf(yb ~ gp(x)), bernoulli()),
  h1 = list(bf(yg ~ gp(x, k = 15)), gaussian()),
  gby = list(bf(yg ~ g + gp(x, by = g)), gaussian()),
  hp = list(bf(yp ~ gp(x, k = 12)), poisson()),
  gs = list(bf(yg ~ gp(x), sigma ~ x2), gaussian()))
out <- list()
for (nm in names(rows)) {
  t0 <- proc.time()[["elapsed"]]
  q <- gf_fit(rows[[nm]][[1]], data = d, family = rows[[nm]][[2]])
  w <- proc.time()[["elapsed"]] - t0
  f <- q$fit
  if (inherits(f, "error")) {
    cat("BIT", nm, "ERROR", conditionMessage(f), "\n"); next
  }
  se <- tryCatch(sqrt(diag(vcov(f, full = TRUE))), error = function(e) NA)
  out[[nm]] <- list(est = f$estimates[c("beta", "betad", "theta", "b")],
                    ll = as.numeric(logLik(f)), se = se,
                    conv = f$opt$convergence, warn = q$warn,
                    gp2 = f$opt$gp_second)
  cat(sprintf("BIT %-6s %s ll %.10f conv %d wall %.2f gp_second %s warn [%s]\n",
              nm, arm, out[[nm]]$ll, out[[nm]]$conv, w,
              paste(format(f$opt$gp_second %||% NA, digits = 4),
                    collapse = ","),
              paste(unique(sub(":.*", "", q$warn)), collapse = ";")))
}
saveRDS(out, sprintf("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/reviews/gpfield-review-out/bit-%s.rds", tag))
