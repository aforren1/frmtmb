# test-se-reparam.R's "nearly unidentified nonlinear body" (seed 2,
# A = 0.05): its stop code, conditions and lost SEs, on each build. The
# CI WARN at :554 was this fit's code-8 warning on ubuntu.
#   CODE8_ARM=before|after Rscript dev/code8-aexp.R
source("dev/code8-common.R")
set.seed(2)
n <- 2000
d <- data.frame(x = stats::runif(n, 1, 2), z = stats::runif(n, 0, 10))
for (k in 1:6) d[[paste0("x", k)]] <- d$x^k
d$y <- sin(3 * d$x) + 0.05 * exp(-0.5 * d$z) + stats::rnorm(n, 0, 0.01)
f <- bf(y ~ b0 + b1 * x1 + b2 * x2 + b3 * x3 + b4 * x4 + b5 * x5 +
          b6 * x6 + A * exp(-k * z),
        b0 + b1 + b2 + b3 + b4 + b5 + b6 + A + k ~ 1, nl = TRUE)
l0 <- stats::lm(y ~ x1 + x2 + x3 + x4 + x5 + x6, data = d)
st <- list(beta = c(unname(stats::coef(l0)), 0.05, 0.5))
r <- rfit(f, data = d, start = st)
cat("code", r$f$opt$convergence, r$f$opt$message, "\n")
cat(gsub(" || ", "\n", cond_txt(r), fixed = TRUE), "\n")
print(stop_info(r))
