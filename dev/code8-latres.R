# test-latres.R's cancelling roots (arma, white noise, seed 2): what each
# build says.
#   CODE8_ARM=before|after Rscript dev/code8-latres.R
source("dev/code8-common.R")
set.seed(2L)
da <- expand.grid(t = 1:6, g = factor(1:10))
da$x <- stats::rnorm(60L)
da$y <- stats::rpois(60L, exp(0.5 + 0.3 * da$x +
                                rep(stats::rnorm(10L, 0, 0.5), each = 6L)))
r <- rfit(bf(y ~ x + arma(t, g, cov = TRUE)), family = poisson(), data = da)
cat("code", r$f$opt$convergence, r$f$opt$message, "\n")
cat(gsub(" || ", "\n", cond_txt(r), fixed = TRUE), "\n")
for (m in r$m) cat("MSG:", m, "\n")
for (w in r$w) cat("WARN:", w, "\n")
print(stop_info(r))
sdr <- quiet(fns$sdr_of(r$f))
print(data.frame(name = fns$outer_par_names(r$f), est = r$f$opt$par,
                 se = sqrt(diag(sdr$cov.fixed))), digits = 4)
