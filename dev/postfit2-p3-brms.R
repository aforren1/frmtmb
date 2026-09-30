# Punch round 3 (P2-M1): mm(g1, g2) + (1 | g1) on draws against brms
# 2.23.0 at the same draws. Two new terms are drawn, and the packages
# consume random numbers in different orders (brms: term outer, draw
# inner; frmtmb: draw outer), so matched seeds do not give equal bands.
# The check is of the law instead: per draw, the new-level offset at
# the middle grid point (spaghetti under re_formula = NULL minus
# spaghetti under NA) divided by its model SD must have variance 1.
#   Rscript dev/postfit2-p3-brms.R [R] > dev/postfit2-log/p3-brms[-R].txt
# Data seed 61 (the review's shape C), draws seed 1, curve seeds 1:R.
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-postfit2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(brms)
  library(frmtmb)
  library(frmtmb.sample)
})
options(mc.cores = 1)
say <- function(...) cat(sprintf(...), "\n", sep = "")
arr <- function(v) array(v, dim = length(v))
N <- 100
R <- if (length(commandArgs(TRUE))) as.integer(commandArgs(TRUE)[1]) else 30
set.seed(61)
n <- 240
d <- data.frame(x = rnorm(n), g1 = factor(sample(1:10, n, TRUE)),
                g2 = factor(sample(1:10, n, TRUE)),
                g3 = factor(sample(1:10, n, TRUE)))
d$w1 <- runif(n, 0.2, 0.8)
d$w2 <- 1 - d$w1
u <- rnorm(10, 0, 1)
v <- rnorm(10, 0, 0.7)
d$y <- rnorm(n, 1 + 0.5 * d$x + d$w1 * u[d$g1] + d$w2 * u[d$g2], 0.5)
d$y3 <- rnorm(n, 1 + 0.5 * d$x + (u[d$g1] + u[d$g2] + u[d$g3]) / 3, 0.5)
d$yp <- rnorm(n, 1 + 0.5 * d$x + 0.5 * (u[d$g1] + u[d$g2]) + v[d$g1], 0.5)
fit <- frm(bf(yp ~ x + (1 | mm(g1, g2)) + (1 | g1)), family = gaussian(),
           data = d)
tpl <- fit$frame$par_template
est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
set.seed(1)
M <- matrix(rep(est, each = N) + rnorm(N * length(est), 0, 0.05), N,
            dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
ds <- structure(list(stanfit = NULL, draws = M, fit = fit),
                class = "frmtmb_draws")
zm <- function(i, grp, s) {
  matrix(M[i, sprintf("r_%s[%d,Intercept]", grp, 1:10)] / s, 1)
}
# theta_1 is the mm() term's, theta_2 (1 | g1)'s; brms numbers g1 first
s_mm <- exp(M[, "theta_1"])
s_g1 <- exp(M[, "theta_2"])
inits <- lapply(seq_len(N), function(i) {
  list(b = arr(M[i, "b_x"]),
       Intercept = M[i, "b_Intercept"] + mean(d$x) * M[i, "b_x"],
       sigma = M[i, "sigma"],
       sd_2 = arr(s_mm[i]), z_2 = zm(i, "mmg1g2", s_mm[i]),
       sd_1 = arr(s_g1[i]), z_1 = zm(i, "g1", s_g1[i]))
})
bo <- suppressMessages(suppressWarnings(
  brm(yp ~ x + (1 | mm(g1, g2)) + (1 | g1), data = d,
      algorithm = "fixed_param", chains = N, iter = 1, warmup = 0,
      init = inits, refresh = 0, seed = 1, silent = 2)))
say("brms sd_mmg1g2 / sd_g1 max rel diff %.3g %.3g",
    max(abs(as_draws_matrix(bo, variable = "sd_mmg1g2__Intercept") /
              s_mm - 1)),
    max(abs(as_draws_matrix(bo, variable = "sd_g1__Intercept") /
              s_g1 - 1)))
spag <- function(o, cond, ...) {
  s <- attr(suppressMessages(suppressWarnings(conditional_effects(
    o, "x", resolution = 3, spaghetti = TRUE, conditions = cond,
    ...)))[[1]], "spaghetti")
  s$estimate__[seq(2, nrow(s), by = 3)]
}
pop_b <- spag(bo, NULL)
pop_f <- spag(ds, list())
say("population spaghetti, brms vs frmtmb: max diff %.3g",
    max(abs(pop_b - pop_f)))
# the model variance of the offset per draw: w_mm^2 * sd_mm^2 summed
# over the distinct new mm levels, plus sd_g1^2
law <- list(
  "nothing set" = list(cb = NULL, cf = list(),
                       v = s_mm^2 + s_g1^2,
                       vb = s_mm^2 + s_g1^2),
  "g1 = 99, g2 unset" = list(cb = data.frame(g1 = "99"),
                             cf = list(g1 = "99"),
                             v = 0.5 * s_mm^2 + s_g1^2,
                             vb = s_mm^2 + s_g1^2))
for (nm in names(law)) {
  L <- law[[nm]]
  zb <- zf <- numeric(0)
  for (r in seq_len(R)) {
    set.seed(r)
    ub <- spag(bo, L$cb, re_formula = NULL, sample_new_levels = "gaussian") -
      pop_b
    uf <- spag(ds, L$cf, re_formula = NULL, seed = r) - pop_f
    zb <- c(zb, ub / sqrt(L$v))
    zf <- c(zf, uf / sqrt(L$v))
  }
  say("%s: %d offsets per package; var(offset / frmtmb's model SD): brms %.4f, frmtmb %.4f (SE about %.4f); brms's own model: %.4f",
      nm, length(zb), var(zb), var(zf), sqrt(2 / length(zb)),
      var(zb * sqrt(L$v / L$vb)))
}
say("done")
