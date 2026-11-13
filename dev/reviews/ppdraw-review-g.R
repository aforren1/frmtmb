# Review G: (1) REML se(sigma = TRUE); (2) message per call; (3) group sd
# near zero over seeds: which interior estimates (lme4 not singular) does
# the scan hold; (4) a quasi-separated logistic slope.
# Usage: Rscript ppdraw-review-g.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
source("C:/Users/adf44/source/r/frmtmb-wt-ppdraw/dev/reviews/ppdraw-review-common.R")
suppressMessages(library(lme4))

cat("\n== G1 REML se(sigma = TRUE)\n")
set.seed(7)
dg <- data.frame(x = stats::rnorm(200), s = stats::runif(200, 0.2, 1))
dg$y <- 1 + dg$x + stats::rnorm(200, 0, dg$s)
fr <- q(frm(bf(y | se(s, sigma = TRUE) ~ x), family = gaussian(), data = dg,
            REML = TRUE))
nd <- data.frame(x = c(0, 0), s = c(0.05, 5))
r <- pz(fr, newdata = nd, ndraws = 400)
cat("held:", length(r$held) > 0, " col sds:",
    format(apply(r$p, 2, stats::sd), digits = 4), "\n")
if (length(r$warns)) cat("warnings:", substr(r$warns, 1, 120), sep = "\n  ")
sim <- as.matrix(simulate(fr, nsim = 400, newdata = nd, seed = 7))
cat("simulate sds:", format(apply(sim, 1, stats::sd), digits = 4), "\n")

cat("\n== G2 messages per call (3269 fit, three predict calls)\n")
set.seed(7)
dz <- data.frame(x = stats::rnorm(200), N = sample(5:15, 200, TRUE))
dz$y <- stats::rbinom(200, dz$N, stats::plogis(0.3 + 0.5 * dz$x))
fz <- q(frm(bf(y | trials(N) ~ x), family = zero_inflated_binomial(),
            data = dz))
k <- 0L
for (i in 1:3) withCallingHandlers(predict(fz, ndraws = 5),
  frmtmb_draw_held = function(m) {
    k <<- k + 1L
    invokeRestart("muffleMessage")
  })
cat("frmtmb_draw_held messages over 3 calls:", k, "\n")
k <- 0L
invisible(tryCatch(withCallingHandlers(predict(fz, ndraws = 5, re_formula = NA),
  frmtmb_draw_held = function(m) {
    k <<- k + 1L
    invokeRestart("muffleMessage")
  }), error = function(e) cat("error:", conditionMessage(e), "\n")))
cat("with re_formula = NA:", k, "\n")

cat("\n== G3 group sd near zero, 6 groups x 20, true sd 0.1, seeds 1..40\n")
cat(sprintf("%4s %8s %8s %10s %10s %8s %8s %5s\n", "seed", "theta",
            "sd", "minus2sd", "plus2sd", "lme4sd", "singular", "held"))
for (seed in 1:40) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(1:6, each = 20)), x = stats::rnorm(120))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(6, 0, 0.1)[d$g] + stats::rnorm(120)
  fit <- q(frm(bf(y ~ x + (1 | g)), data = d))
  fs <- flat_sides(fit)
  if (is.null(fs)) {
    cat(sprintf("%4d no usable covariance\n", seed))
    next
  }
  if (seed == 1) print(rownames(fs))
  j <- which(!grepl("^\\(Intercept\\)$|^x$|^sigma", rownames(fs)))[1]
  lm4 <- suppressMessages(lmer(y ~ x + (1 | g), data = d, REML = FALSE))
  r <- pz(fit, ndraws = 20)
  m <- min(abs(fs[j, c("minus2sd", "plus2sd")]))
  if (isTRUE(m < 0.05) || length(r$held)) {
    cat(sprintf("%4d %8.3f %8.3g %10.3g %10.3g %8.4f %8s %5s\n", seed,
                fs[j, "est"], fs[j, "sd"], fs[j, "minus2sd"], fs[j, "plus2sd"],
                attr(VarCorr(lm4)$g, "stddev"), isSingular(lm4),
                length(r$held) > 0))
  }
}

cat("\n== G4 quasi-separated logistic slope\n")
x <- seq(-3, 3, length.out = 80)
y <- as.integer(x > 0)
y[x > -0.2 & x < 0] <- 1L
y[x > 0 & x < 0.2] <- 0L
dl <- data.frame(x = x, y = y)
fl <- q(frm(bf(y ~ x), family = bernoulli(), data = dl))
print(signif(flat_sides(fl), 4))
r <- pz(fl, ndraws = 400, newdata = data.frame(x = c(-0.3, 0, 0.3)))
cat("held:", length(r$held) > 0, " predict means:",
    format(colMeans(r$p), digits = 4), " fitted:",
    format(fitted(fl, newdata = data.frame(x = c(-0.3, 0, 0.3)))[, "Estimate"],
           digits = 4), "\n")
