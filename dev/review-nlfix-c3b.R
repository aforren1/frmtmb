# Item 2(b), 3: flagged fits that restarts cannot help, large enough to feel. Seed 6403.
source("dev/review-nlfix-common.R")
q <- function(expr) suppressWarnings(suppressMessages(expr))
set.seed(6403)
n <- 2e4; ng <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(sample(ng, n, TRUE)))
d$y <- 2 + 0.5 * d$x + rnorm(ng, 0, 0.3)[as.integer(d$g)] + rnorm(n, 0, 0.5)
d$yo <- cut(0.8 * d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = 1:4, ordered_result = TRUE)
tm <- function(label, default, zero) {
  a <- system.time(ra <- default())[["elapsed"]]; a2 <- system.time(default())[["elapsed"]]
  b <- system.time(rb <- zero())[["elapsed"]]; b2 <- system.time(zero())[["elapsed"]]
  cat(sprintf("%-32s default %.2f s (ll %.4f conv %s), one fit %.2f s (ll %.4f conv %s), ratio %.1f\n", label, min(a, a2), as.numeric(logLik(ra)), ra$opt$convergence, min(b, b2), as.numeric(logLik(rb)), rb$opt$convergence, min(a, a2) / min(b, b2)))
}
z <- function(nm) list(beta = setNames(rep(0, length(nm)), nm))
f1 <- bf(y ~ a + b + c0 * x, a ~ 1 + (1 | g), b ~ 1, c0 ~ 1, nl = TRUE)
nm1 <- c("a_(Intercept)", "b_(Intercept)", "c0_(Intercept)")
tm("RE, ridge a + b", function() q(frm(f1, data = d)), function() q(frm(f1, data = d, start = z(nm1))))
f2 <- bf(yo ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)
nm2 <- c("b0_(Intercept)", "b1_(Intercept)")
tm("cumulative, ridge b0 + b1 x", function() q(frm(f2, family = cumulative(), data = d)), function() q(frm(f2, family = cumulative(), data = d, start = z(nm2))))
f3 <- bf(y ~ a * exp(b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE)
ct <- frmtmb_control(optCtrl = list(iter.max = 3L))
nm3 <- c("a_(Intercept)", "b_(Intercept)")
tm("RE, iter.max = 3 (not converged)", function() q(frm(f3, data = d, control = ct)), function() q(frm(f3, data = d, control = ct, start = z(nm3))))
