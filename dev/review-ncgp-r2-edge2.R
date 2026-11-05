# Fix round 1, item 1: gaussian spellings the gate turns away, on 2 sin(3x) + x. Default against the best of a
# grid of starts on the same tape. Args: arm. Seed 58001.
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs(arm)
set.seed(58001); x <- sort(runif(60, 0, 10)); z <- rnorm(60); zm <- z; zm[seq(5, 60, 12)] <- NA
d <- data.frame(x, z, zm, zo = z + rnorm(60, sd = 0.1), s = 0.1, y = 2 * sin(3 * x) + x + 0.3 * z + 0.1 * rnorm(60))
d$o <- ordered(sample(1:4, 60, TRUE))
sp <- list(
  nointercept = list(y ~ 0 + gp(x)), start_beta = list(y ~ gp(x), start = list(beta = 5)), start_bd = list(y ~ gp(x), start = list(betad = log(0.5))),
  mi = list(bf(y ~ mi(zm) + gp(x)) + bf(zm | mi() ~ 1) + set_rescor(FALSE)), me = list(y ~ me(zo, s) + gp(x)), mo = list(y ~ mo(o) + gp(x)),
  cs = list(y ~ gp(x) + cs(0 + factor(as.integer(o)) | rep(1:6, 10))), sigma_gp = list(bf(y ~ gp(x), sigma ~ gp(z))), nl = list(bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE)),
  mix = list(bf(y ~ gp(x)), family = mixture(gaussian(), gaussian()))
)
for (nm in names(sp)) {
  s <- sp[[nm]]
  t0 <- proc.time()[[3]]; r <- do.call(rev_fit, c(s[1], list(data = d), s[-1])); el <- proc.time()[[3]] - t0
  if (inherits(r$fit, "error")) { cat(sprintf("RE3 %s %-11s ERROR %s\n", arm, nm, substr(conditionMessage(r$fit), 1, 110))); next }
  ll0 <- as.numeric(logLik(r$fit)); best <- ll0; th0 <- r$fit$estimates$theta
  for (rho0 in c(0.15, 0.5, 1.5, 5)) for (sd0 in c(0.5, 3, 15)) {
    th <- th0
    for (b in r$fit$frame$re_blocks) if (identical(b$covstruct, "gp")) { th[b$theta_idx[1]] <- log(sd0); th[b$theta_idx[-1]] <- log(rho0) }
    a <- s[-1]; a$start <- c(a$start, list(theta = th))
    r2 <- do.call(rev_fit, c(s[1], list(data = d), a))
    if (!inherits(r2$fit, "error")) { v <- as.numeric(logLik(r2$fit)); if (is.finite(v) && v > best) best <- v }
  }
  cat(sprintf("RE3 %s %-11s default ll %.4f conv %d w %d | best %.4f | gap %.2e %.1fs\n", arm, nm, ll0, r$fit$opt$convergence, length(r$warn), best, ll0 - best, el))
  flush.console()
}
cat("DONE\n")
