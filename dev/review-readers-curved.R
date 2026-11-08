# Review of lane readers, short check of M3: conditional_effects() on two
# more bodies whose nlpars read me(), against fitted() and Monte Carlo
# over the draw at the estimates (seed 81, 1e6 draws).
# A: a steep sigmoid, y ~ 10 / (1 + exp(-4 a)), a ~ me(xo, sx); data seed
#    8301. B: one draw read by two nlpars, y ~ a * exp(b), a ~ me(xo, sx),
#    b ~ me(xo, sx); data seed 8302.
source("dev/review-readers-setup.R")
dm <- function(f, nd, dp) as.numeric(frm_linpred(f, newdata = nd, dpar = dp))
slope <- function(f, nd, dp, v) {
  nd1 <- nd
  nd1[[v]] <- nd1[[v]] + 1
  dm(f, nd1, dp) - dm(f, nd, dp)
}
n <- 600
report <- function(f, body, sds, x0) {
  set.seed(81)
  z <- rnorm(1e6)
  for (s in sds) {
    nd <- data.frame(xo = x0, sx = s)
    r <- with_warn(fitted(f, newdata = nd))
    ce <- with_warn(conditional_effects(f, effects = "xo",
                                        int_conditions = list(xo = x0),
                                        conditions = data.frame(sx = s)))
    c1 <- ce$value[[1]]
    mc <- body(nd, z)
    cat(sprintf(paste("sx %.2f: fitted %.4f (%.4f) | ce est %.4f se %.4f",
                      "band %.4f to %.4f | MC %.4f (%.4f) q %.4f to %.4f |",
                      "warn %s\n"),
                s, r$value[, "Estimate"], r$value[, "Est.Error"],
                c1$estimate__, c1$se__, c1$lower__, c1$upper__, mean(mc),
                sd(mc), quantile(mc, 0.025), quantile(mc, 0.975),
                paste(unique(c(wclass(r$warnings), wclass(ce$warnings))),
                      collapse = ",")))
  }
}

cat("== A. y ~ 10 / (1 + exp(-4 a)), a ~ me(xo, sx)\n")
set.seed(8301)
d <- data.frame(l = rnorm(n), sx = 0.3)
d$xo <- d$l + rnorm(n, 0, 0.3)
d$y <- 10 / (1 + exp(-4 * (0.2 + 0.5 * d$l))) + rnorm(n, 0, 0.5)
fa <- frm(bf(y ~ 10 / (1 + exp(-4 * a)), a ~ me(xo, sx), nl = TRUE),
          family = gaussian(), data = d)
ba <- slope(fa, data.frame(xo = 0, sx = 1), "a", "xo")
cat("slope of a:", ba, "\n")
report(fa, function(nd, z) {
  10 / (1 + exp(-4 * (dm(fa, nd, "a") + ba * nd$sx * z)))
}, c(0.3, 1, 2), 0)

cat("\n== B. y ~ a * exp(b), a and b read one me() call\n")
set.seed(8302)
d <- data.frame(l = rnorm(n), sx = 0.3)
d$xo <- d$l + rnorm(n, 0, 0.3)
d$y <- (2 + 0.5 * d$l) * exp(0.3 + 0.4 * d$l) + rnorm(n, 0, 0.3)
fb <- frm(bf(y ~ a * exp(b), a ~ me(xo, sx), b ~ me(xo, sx), nl = TRUE),
          family = gaussian(), data = d)
b_a <- slope(fb, data.frame(xo = 0, sx = 1), "a", "xo")
b_b <- slope(fb, data.frame(xo = 0, sx = 1), "b", "xo")
cat("slopes a, b:", b_a, b_b, "\n")
report(fb, function(nd, z) {
  w <- nd$sx * z
  (dm(fb, nd, "a") + b_a * w) * exp(dm(fb, nd, "b") + b_b * w)
}, c(0.3, 1, 2), 0.5)
