# Lane ceplot: fitted(sample_new_levels = "old_levels") on a fit.
# 1. The answer at an unseen level IS the fitted answer at one seen
#    level (estimate and Est.Error, compared with identical()).
# 2. The level is chosen as brms chooses it: at the same seed, brms
#    2.23.0's fitted() on draws of this model picks the same level
#    (brms at algorithm = "fixed_param", initialized at draws around
#    the fit, as dev/postfit2-p2-brms.R does).
# 3. predict() at the same seed picks the same level (the shared code):
#    read from the mean of its draws, and from its own choice
#    (predict_new_level_spec()), since two seen levels can lie closer
#    together than the Monte Carlo error of the mean.
#   Rscript dev/ceplot-oldlevels.R > dev/ceplot-log/oldlevels.txt 2>&1
# Seeds: data 3, draws 1, calls 1 to 12.
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(brms)
  library(frmtmb)
})
cat("frmtmb from", find.package("frmtmb"), "\n")
say <- function(...) cat(sprintf(...), "\n", sep = "")
set.seed(3)
dg <- data.frame(x = rnorm(120), g = factor(rep(1:12, 10)))
dg$y <- rnorm(120, 1 + 0.5 * dg$x + rnorm(12, 0, 2)[dg$g])
fg <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dg)
nd <- data.frame(x = c(0, 1), g = factor(c("new", "new")))
at <- lapply(1:12, function(k) {
  fitted(fg, newdata = data.frame(x = c(0, 1), g = factor(k, levels = 1:12)))
})
pick_fit <- integer(12)
for (s in 1:12) {
  set.seed(s)
  f <- fitted(fg, newdata = nd, allow_new_levels = TRUE,
              sample_new_levels = "old_levels")
  hit <- which(vapply(at, function(a) {
    identical(unname(a[, 1:2]), unname(f[, 1:2]))
  }, NA))
  pick_fit[s] <- if (length(hit) == 1L) hit else NA_integer_
}
say("fitted: seeds 1..12 read level %s (identical() to fitted at it)",
    paste(pick_fit, collapse = " "))
gauss <- fitted(fg, newdata = nd, allow_new_levels = TRUE)
say("Est.Error: gaussian %.4f %.4f | old_levels at seed 1 %.4f %.4f",
    gauss[1, 2], gauss[2, 2], at[[pick_fit[1]]][1, 2],
    at[[pick_fit[1]]][2, 2])
cand <- fixef(fg)["Intercept", "Estimate"] + ranef(fg)$g[, "Intercept"]
pick_pred <- integer(12)
for (s in 1:12) {
  set.seed(s)
  d <- predict(fg, newdata = nd[1, ], allow_new_levels = TRUE,
               sample_new_levels = "old_levels", propagate_error = FALSE,
               ndraws = 20000, summary = FALSE)
  pick_pred[s] <- which.min(abs(mean(d) - cand))
}
say("predict: seeds 1..12 read level %s", paste(pick_pred, collapse = " "))
say("seen levels (fixef + ranef): %s", paste(sprintf("%.3f", cand), collapse = " "))
say("predict and fitted pick the same level at %d of 12 seeds", sum(pick_pred == pick_fit))
pick_spec <- integer(12)
for (s in 1:12) {
  set.seed(s)
  sp <- frmtmb:::predict_new_level_spec(fg, fg$spec$responses[[1L]], nd,
                                        NULL, TRUE, "old_levels")
  pick_spec[s] <- attr(sp, "old_pick")[[1L]]
}
say("predict's own choice (predict_new_level_spec): seeds 1..12 level %s",
    paste(pick_spec, collapse = " "))
say("predict's choice and fitted's agree at %d of 12 seeds",
    sum(pick_spec == pick_fit))

## brms at draws around the fit
N <- 50
set.seed(1)
inits <- lapply(seq_len(N), function(i) {
  sd1 <- VarCorr(fg)$g$sd[1, 1] * exp(rnorm(1, 0, 0.05))
  r <- ranef(fg)$g[, "Intercept"] + rnorm(12, 0, 0.05)
  list(b = array(fixef(fg)["x", "Estimate"] + rnorm(1, 0, 0.05), 1),
       Intercept = fixef(fg)["Intercept", "Estimate"] +
         mean(dg$x) * fixef(fg)["x", "Estimate"] + rnorm(1, 0, 0.05),
       sigma = sigma(fg), sd_1 = array(sd1, 1), z_1 = matrix(r / sd1, 1))
})
b <- suppressMessages(suppressWarnings(
  brm(y ~ x + (1 | g), data = dg, algorithm = "fixed_param", chains = N,
      iter = 1, warmup = 0, init = inits, refresh = 0, seed = 1,
      silent = 2)))
bat <- lapply(1:12, function(k) {
  fitted(b, newdata = data.frame(x = c(0, 1), g = factor(k, levels = 1:12)),
         summary = FALSE)
})
pick_brms <- integer(12)
for (s in 1:12) {
  set.seed(s)
  f <- fitted(b, newdata = nd, allow_new_levels = TRUE,
              sample_new_levels = "old_levels", summary = FALSE)
  hit <- which(vapply(bat, function(a) isTRUE(all.equal(a, f)), NA))
  pick_brms[s] <- if (length(hit) == 1L) hit else NA_integer_
}
say("brms: seeds 1..12 read level %s", paste(pick_brms, collapse = " "))
say("brms and frmtmb fitted pick the same level at %d of 12 seeds",
    sum(pick_brms == pick_fit, na.rm = TRUE))
say("done")
