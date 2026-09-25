source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-helpers.R")
# Does the default engage on ordinary data (slope column sd in 0.005 to
# 0.045) and, when it does, does it reproduce autoscale = FALSE? Also
# the singular case (true slope sd 0). Seeds 101..104.
ctlF <- frmtmb_control(autoscale = FALSE)
res <- list()
for (fam in c("gaussian", "poisson", "bernoulli"))
for (sx in c(0.005, 0.02, 0.045)) for (ssd in c(0, 0.3)) for (seed in 101:104) {
  set.seed(seed)
  ng <- 25; per <- 12; n <- ng * per
  g <- factor(rep(seq_len(ng), each = per))
  xs <- runif(n, 0, 1); x <- xs * sx / sd(xs)   # a proportion-like column
  u0 <- rnorm(ng, 0, 0.6); u1 <- rnorm(ng, 0, ssd)
  eta <- 0.2 + 0.5 * xs + u0[g] + u1[g] * xs
  y <- switch(fam, gaussian = eta + rnorm(n),
              poisson = rpois(n, exp(eta)),
              bernoulli = rbinom(n, 1, plogis(eta)))
  d <- data.frame(y, x, g)
  famo <- switch(fam, gaussian = gaussian(), poisson = poisson(),
                 bernoulli = bernoulli())
  t1 <- proc.time()[3]
  a <- fitw(y ~ x + (1 + x | g), data = d, family = famo)
  t2 <- proc.time()[3]
  b <- fitw(y ~ x + (1 + x | g), data = d, family = famo, control = ctlF)
  t3 <- proc.time()[3]
  la <- if (inherits(a$fit, "error")) NA else as.numeric(logLik(a$fit))
  lb <- if (inherits(b$fit, "error")) NA else as.numeric(logLik(b$fit))
  vca <- tryCatch(vc_num(a$fit), error = function(e) NA)
  vcb <- tryCatch(vc_num(b$fit), error = function(e) NA)
  res[[length(res) + 1L]] <- data.frame(
    fam, sx, ssd, seed, engaged = !is.null(a$tpl), dLL = la - lb,
    llF = lb, codeD = if (is.na(la)) NA else a$fit$opt$convergence,
    codeF = if (is.na(lb)) NA else b$fit$opt$convergence,
    vc_rel = rel(vca, vcb), fe_rel = if (is.na(la) || is.na(lb)) NA else
      rel(fixef(a$fit), fixef(b$fit)),
    nwD = length(a$warn), nwF = length(b$warn),
    warnsame = identical(sort(unique(sub("[(].*", "", a$warn))),
                         sort(unique(sub("[(].*", "", b$warn)))),
    tD = t2 - t1, tF = t3 - t2)
}
r <- do.call(rbind, res)
options(width = 200)
print(r, digits = 4)
cat("\nfits:", nrow(r), " engaged:", sum(r$engaged),
    " default lower by >1e-6:", sum(r$dLL < -1e-6, na.rm = TRUE),
    " higher by >1e-6:", sum(r$dLL > 1e-6, na.rm = TRUE),
    " max |dLL|:", format(max(abs(r$dLL), na.rm = TRUE), digits = 3),
    " warning sets differ:", sum(!r$warnsame),
    " code differs:", sum(r$codeD != r$codeF, na.rm = TRUE), "\n")
cat("time ratio default/F, engaged: median",
    format(median(r$tD[r$engaged] / r$tF[r$engaged]), digits = 3), "\n")
saveRDS(r, "C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-m1-battery.rds")
