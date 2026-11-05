# The ordinal mean display against brms, body and linear; the patched density's cost.
source("dev/review-nlfix-common.R")
suppressPackageStartupMessages(library(testthat))
ns <- asNamespace("frmtmb")
E <- new.env(parent = ns); E$`%||%` <- function(a, b) if (is.null(a)) b else a
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) sys.source(f, envir = E)
set.seed(9202); n <- 300
d <- data.frame(x = round(rnorm(n), 3)); d$y <- ordered(cut(1.5 * exp(0.6 * d$x) + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
r <- suppressWarnings(suppressMessages(frm(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), family = cumulative(), data = d)))
l <- frm(y ~ x, family = cumulative(), data = d)
fr <- suppressWarnings(E$brms_fixed_fit(brms::bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), brms::cumulative(), d, r))
fl <- suppressWarnings(E$brms_fixed_fit(brms::bf(y ~ x), brms::cumulative(), d, l))
a <- conditional_effects(r, effects = "x", resolution = 5)[[1]]; b <- suppressWarnings(brms::conditional_effects(fr, effects = "x", resolution = 5)[[1]])
print(a[, c("x", "estimate__", "se__")]); print(b[, c("x", "estimate__", "se__")])
a2 <- conditional_effects(l, effects = "x", resolution = 5)[[1]]; b2 <- suppressWarnings(brms::conditional_effects(fl, effects = "x", resolution = 5)[[1]])
print(a2[, c("x", "estimate__")]); print(b2[, c("x", "estimate__")])
cat("frm body, what is it: sum k p_k with k = 1..4:\n")
p <- fitted(r, newdata = data.frame(x = a$x))[, "Estimate", ]
print(round(p %*% 1:4, 5))
cat("\ncost of the patched density, cumulative n = 3e4 row lpdf in a fit\n")
set.seed(21); nn <- 3e4; dd <- data.frame(x = rnorm(nn), ntr = 10L); p <- plogis(0.3 + 0.5 * dd$x); dd$y <- rbinom(nn, 10, rbeta(nn, p * 8, (1 - p) * 8))
bb2 <- frmtmb_family("bb2", dpars = c("mu", "phi"), links = list(mu = "logit", phi = "log"), type = "discrete", required_aterms = "trials",
  lpdf = function(y, dpars, aterms = NULL, extra = NULL) {
    size <- ns$aterm_trials(aterms, "bb2"); a <- dpars$mu * dpars$phi; b <- (1 - dpars$mu) * dpars$phi
    lgamma(size + 1) - lgamma(y + 1) - lgamma(size - y + 1) + ns$lgamma_shift_diff(a, y) + ns$lgamma_shift_diff(b, size - y) - ns$lgamma_shift_diff(dpars$phi, size)
  }, init_dpars = list(mu = function(y, aterms) 0.5, phi = function(y, aterms) 5))
t1 <- Inf; for (i in 1:3) t1 <- min(t1, system.time(f1 <- frm(bf(y | trials(ntr) ~ x), family = beta_binomial(), data = dd))[["elapsed"]])
t2 <- Inf; for (i in 1:3) t2 <- min(t2, system.time(f2 <- frm(bf(y | trials(ntr) ~ x), family = bb2, data = dd))[["elapsed"]])
cat(sprintf("current %.3f s ll %.6f | patched %.3f s ll %.6f | ratio %.2f\n", t1, as.numeric(logLik(f1)), t2, as.numeric(logLik(f2)), t2 / t1))
d3 <- dd[1:2000, ]; d3$y <- pmin(d3$y, 10); d3$yk <- rbinom(2000, 200, plogis(0.3 + 0.5 * d3$x)); d3$nk <- 200L; d3$yz <- ifelse(runif(2000) < .1, 0L, d3$yk)
cat("trials = 200: current", format(as.numeric(logLik(suppressWarnings(frm(bf(yk | trials(nk) ~ x), family = beta_binomial(), data = d3)))), digits = 10), " patched", format(as.numeric(logLik(suppressWarnings(frm(bf(yk | trials(nk) ~ x), family = bb2, data = d3)))), digits = 10), " binomial", format(as.numeric(logLik(glm(cbind(yk, nk - yk) ~ x, binomial, d3))), digits = 10), "\n")
cat("\nzero_inflated_beta_binomial at the same large phi: objective on doubles\n")
d9 <- { set.seed(9); n <- 300; x <- data.frame(x = round(rnorm(n), 3), ntr = 10L); x$y <- rbinom(n, 10, plogis(0.3 + 0.5 * x$x)); x }
fz <- suppressWarnings(frm(bf(y | trials(ntr) ~ x), family = zero_inflated_beta_binomial(), data = d9))
print(fz$estimates$betad); tp <- fz$estimates; nll <- ns$build_objective(fz$frame)
for (lp in c(20, 30, 39, 45)) { tp$betad[1] <- lp; cat(lp, -as.numeric(nll(tp)), "\n") }
