# brms 2.23.0's posterior SD of the quantities dev/mfxse-rate.R measures, on a
# subset of its designs (seed 7501), beside the marginaleffects SE of each arm
# read from dev/mfxse-rate-out/{before,after}-<family>.tsv. brms's default
# priors; 2 chains of 2000 iterations, run one after the other.
#   Rscript dev/mfxse-brms.R
source("C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/mfxse-lib.R")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
O <- file.path(mfx_root, "dev/mfxse-rate-out/")
subset <- data.frame(
  family = c("poisson", "poisson", "gaussian", "bernoulli", "lognormal"),
  J = c(30L, 30L, 30L, 30L, 30L), m = c(10L, 10L, 10L, 40L, 10L),
  sdg = 0.8, slope = c(FALSE, TRUE, FALSE, FALSE, FALSE), seed = 7501L)
gen <- function(famname, eta, n) {
  switch(famname,
         gaussian = eta + rnorm(n, 0, 0.8),
         poisson = rpois(n, exp(eta)),
         bernoulli = rbinom(n, 1, plogis(eta)),
         lognormal = exp(eta + rnorm(n, 0, 0.5)))
}
bfam <- list(gaussian = gaussian(), poisson = poisson(),
             bernoulli = bernoulli(), lognormal = lognormal())
out <- list()
for (i in seq_len(nrow(subset))) {
  s <- subset[i, ]
  set.seed(s$seed)
  n <- s$J * s$m
  d <- data.frame(g = factor(rep(sprintf("g%02d", seq_len(s$J)), each = s$m)),
                  x = rnorm(n), f = factor(sample(c("a", "b"), n, TRUE)))
  gw <- rnorm(s$J)
  d$w <- gw[d$g]
  u0 <- rnorm(s$J, 0, s$sdg)[d$g]
  u1 <- if (s$slope) rnorm(s$J, 0, s$sdg / 2)[d$g] else 0
  eta <- 0.4 + (0.5 + u1) * d$x + 0.4 * (d$f == "b") + 0.3 * d$w + u0
  d$y <- gen(s$family, eta, n)
  rhs <- if (s$slope) "x + f + w + (1 + x | g)" else "x + f + w + (1 | g)"
  fml <- stats::as.formula(paste("y ~", rhs))
  t0 <- proc.time()[["elapsed"]]
  bfit <- quiet(brm(fml, family = bfam[[s$family]], data = d, chains = 2,
                    iter = 2000, cores = 1, seed = s$seed, refresh = 0,
                    backend = "rstan"))
  secs <- proc.time()[["elapsed"]] - t0
  ndp <- d[seq(1, n, length.out = 12), ]
  h <- 1e-4
  ep <- function(nd) posterior_epred(bfit, newdata = nd)
  dh <- transform(d, x = x + h)
  dl <- transform(d, x = x - h)
  da <- transform(d, f = factor("a", levels = levels(d$f)))
  db <- transform(d, f = factor("b", levels = levels(d$f)))
  dw1 <- transform(d, w = w + 1)
  sd_ <- function(m) stats::sd(m)
  bsd <- c(apply(ep(ndp), 2, sd),
           avg_pred = sd_(rowMeans(ep(d))),
           slope_x = sd_(rowMeans((ep(dh) - ep(dl)) / (2 * h))),
           cmp_f = sd_(rowMeans(ep(db) - ep(da))),
           cmp_w = sd_(rowMeans(ep(dw1) - ep(d))))
  rh <- max(rhat(bfit), na.rm = TRUE)
  arms <- lapply(c("before", "after"), function(a) {
    r <- read.delim(file.path(O, paste0(a, "-", s$family, ".tsv")))
    r <- r[r$J == s$J & r$m == s$m & r$sdg == s$sdg & r$slope == s$slope &
             r$seed == s$seed, ]
    r
  })
  qty <- arms[[2]]$qty
  res <- data.frame(family = s$family, J = s$J, m = s$m, slope = s$slope,
                    seed = s$seed, qty = qty, brms_sd = unname(bsd),
                    before = arms[[1]]$me_se, after = arms[[2]]$me_se,
                    draws_ref = arms[[2]]$ref, rhat = rh, brms_secs = secs)
  out[[i]] <- res
  cat(sprintf("%s J%d m%d slope %s: max rhat %.3f, %.0f s\n", s$family, s$J,
              s$m, s$slope, rh, secs))
  agg <- res[res$qty != "pred", ]
  print(agg[, c("qty", "brms_sd", "before", "after", "draws_ref")],
        row.names = FALSE, digits = 4)
  pr <- res[res$qty == "pred", ]
  cat(sprintf("  levels: after/brms median %.3f range %.3f to %.3f; before/brms median %.3f range %.3f to %.3f\n",
              median(pr$after / pr$brms_sd), min(pr$after / pr$brms_sd),
              max(pr$after / pr$brms_sd), median(pr$before / pr$brms_sd),
              min(pr$before / pr$brms_sd), max(pr$before / pr$brms_sd)))
}
all <- do.call(rbind, out)
write.table(all, file.path(O, "brms.tsv"), sep = "\t", row.names = FALSE,
            quote = FALSE)
cat("\nratio to brms's posterior SD, all quantities but levels:\n")
a <- all[all$qty != "pred", ]
cat(sprintf("after: median %.3f range %.3f to %.3f; before: median %.3f range %.3f to %.3f\n",
            median(a$after / a$brms_sd), min(a$after / a$brms_sd),
            max(a$after / a$brms_sd), median(a$before / a$brms_sd),
            min(a$before / a$brms_sd), max(a$before / a$brms_sd)))
cat("DONE\n")
