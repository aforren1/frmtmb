# conditional_effects(re_formula = NULL, conditions = <an observed
# level>) on draws against brms at the same five draws. The draws are
# the hand-built ones of dev/postfit2-celevel.R (lane sampfix's
# construction): the ML estimate plus N(0, 0.05) noise, seed 1; data
# seed 9. brms is sampled with algorithm = "fixed_param", one chain per
# draw, initialized at that draw.
#   Rscript dev/postfit2-celevel-brms.R > dev/postfit2-log/celevel-brms.txt
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-postfit2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(brms)
  library(frmtmb)
  library(frmtmb.sample)
})
say <- function(...) cat(sprintf(...), "\n", sep = "")
set.seed(9)
dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
lab <- frmtmb::brms_par_labels(fit)
tpl <- fit$frame$par_template
est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
set.seed(1)
M <- matrix(rep(est, each = 5) + rnorm(5 * length(est), 0, 0.05), 5,
            dimnames = list(NULL, lab))
M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
full <- structure(list(stanfit = NULL, draws = M, fit = fit),
                  class = "frmtmb_draws")
mx <- mean(dd$x)
# the sd only scales z_1 back to r, so a wrong transform would leave the
# curves unchanged; the r_ comparison below is what carries the check
say("exp(theta) %.6f; the VarCorr() sd:", exp(fit$estimates$theta[1]))
print(VarCorr(fit)$g$sd)
inits <- lapply(seq_len(nrow(M)), function(i) {
  sdg <- exp(M[i, "theta_1"])   # a (1 | g) block holds log sd
  r <- M[i, sprintf("r_g[%d,Intercept]", 1:6)]
  list(b = array(M[i, "b_x"], 1),
       Intercept = M[i, "b_Intercept"] + mx * M[i, "b_x"],
       sd_1 = array(sdg, 1), z_1 = matrix(r / sdg, 1, 6),
       sigma = M[i, "sigma"])
})
bK <- suppressMessages(suppressWarnings(
  brm(y ~ x + (1 | g), data = dd, algorithm = "fixed_param",
      chains = length(inits), iter = 1, warmup = 0, init = inits,
      refresh = 0, seed = 1, silent = 2)))
say("brms draws %d; r_g[2,Intercept] brms vs frmtmb max diff %.3g",
    ndraws(bK), max(abs(as_draws_matrix(bK, variable = "r_g[2,Intercept]") -
                          M[, "r_g[2,Intercept]"])))
for (lev in 1:3) {
  cnd <- data.frame(g = factor(lev, levels = 1:6))
  a <- conditional_effects(bK, "x", re_formula = NULL, conditions = cnd,
                           resolution = 7, spaghetti = TRUE)[[1]]
  b <- conditional_effects(full, "x", re_formula = NULL, conditions = cnd,
                           resolution = 7, spaghetti = TRUE)[[1]]
  sa <- attr(a, "spaghetti")
  sb <- attr(b, "spaghetti")
  say("g = %d: estimate__ max |brms - frmtmb| %.3g, se__ %.3g, lower__ %.3g, upper__ %.3g; spaghetti rows %d / %d, estimate__ %.3g",
      lev, max(abs(a$estimate__ - b$estimate__)),
      max(abs(a$se__ - b$se__)), max(abs(a$lower__ - b$lower__)),
      max(abs(a$upper__ - b$upper__)), nrow(sa), nrow(sb),
      max(abs(sa$estimate__ - sb$estimate__)))
}
