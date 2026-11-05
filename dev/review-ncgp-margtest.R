# Item 8: what the closed form needs: by =, several terms, (1 | g), se(), a covariate. The
# marginal tape from the frame at frm()'s optimum against the z tape's objective; the warm
# start against the default; seconds. A low-noise and a rough design per spec. Seed 6100 + i.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
source(file.path(REV_WT, "dev/review-ncgp-margin.R"))
ns <- asNamespace("frmtmb")
specs <- list(
  plain = y ~ gp(x), cov = y ~ z + gp(x), byf = y ~ gp(x, by = f), byw = y ~ gp(x, by = w), cmcF = y ~ gp(x, by = f, cmc = FALSE),
  two = y ~ gp(x) + gp(z), re = y ~ z + gp(x) + (1 | g), slope = y ~ gp(x) + (z | g), se = y | se(s, sigma = TRUE) ~ gp(x), rep = y ~ gp(xr), grF = y ~ gp(xr, gr = FALSE), md = y ~ gp(x, z),
  aniso = y ~ gp(x, z, iso = FALSE), smooth = y ~ gp(x) + s(z)
)
for (sg in c(0.1, 1e-3)) for (i in seq_along(specs)) {
  set.seed(6100 + i)
  n <- 150
  d <- data.frame(x = runif(n, 0, 8), z = runif(n, 0, 8), g = factor(sample(letters[1:8], n, TRUE)), f = factor(sample(c("u", "v", "w"), n, TRUE)),
                  w = runif(n, 0.5, 2), s = 0.05, xr = rep(sort(runif(30, 0, 8)), each = 5))
  d$y <- 2 * sin(3 * d$x) + d$x + 0.5 * sin(d$z) + sin(3 * d$xr) + sg * rnorm(n) + 0.2 * rnorm(8)[as.integer(d$g)] * (names(specs)[i] %in% c("re", "slope"))
  t0 <- proc.time()[[3]]
  r <- rev_fit(specs[[i]], d)
  t1 <- proc.time()[[3]] - t0
  if (inherits(r$fit, "error")) { cat("RMA", names(specs)[i], "FIT ERROR", substr(conditionMessage(r$fit), 1, 100), "\n"); next }
  f <- r$fit
  if (!marg_applicable(f$frame)) { cat(sprintf("RMA %-7s s=%g not applicable: blocks %s\n", names(specs)[i], sg, paste(vapply(f$frame$re_blocks, function(k) k$covstruct, ""), collapse = ","))); next }
  ms <- tryCatch(marg_start(f$frame), error = function(e) e)
  if (inherits(ms, "error")) { cat("RMA", names(specs)[i], "BUILD ERROR", conditionMessage(ms), "\n"); next }
  tp <- f$estimates[c("beta", "betad", "theta")]
  mv <- tryCatch(ms$nll(tp), error = function(e) NA)
  t0 <- proc.time()[[3]]
  mo <- tryCatch(marg_opt(f$frame), error = function(e) e)
  tm <- proc.time()[[3]] - t0
  if (inherits(mo, "error")) { cat("RMA", names(specs)[i], "OPT ERROR", conditionMessage(mo), "\n"); next }
  t0 <- proc.time()[[3]]
  r2 <- rev_fit(specs[[i]], d, start = mo$start)
  t2 <- proc.time()[[3]] - t0
  f2 <- r2$fit
  cat(sprintf("RMA %-7s s=%g marginal-tape at frm's point: diff %.2e | default ll %.4f conv %d w %d %.2fs | marginal ll %.4f %.2fs, then tape ll %.4f conv %d w %d %.2fs | gain over default %.3f\n",
              names(specs)[i], sg, -mv + f$opt$objective * -1 + 0 * mv, as.numeric(logLik(f)), f$opt$convergence, length(r$warn), t1, mo$ll, tm,
              if (inherits(f2, "error")) NA else as.numeric(logLik(f2)), if (inherits(f2, "error")) NA else f2$opt$convergence, length(r2$warn), t2,
              if (inherits(f2, "error")) NA else as.numeric(logLik(f2)) - as.numeric(logLik(f))))
  flush.console()
}
cat("DONE\n")
