# Reviewer, punch round 1: is frmtmb's new-level draw the same law as
# brms's sample_new_levels = "gaussian"? Same 100 draws in both (the
# lane's construction: ML estimate plus N(0, 0.05), seed 1; data seed
# 9; brms fixed_param, one chain per draw). Per repeat, the spaghetti
# of re_formula = NULL (a new group) minus the spaghetti of re_formula
# = NA is that draw's new-level effect u_i; z_i = u_i / sd_i. The same
# law gives var(z) = 1 in both. R repeats with seeds 1..R.
#   Rscript dev/postfit2-rev-p1-law.R [R]
args <- commandArgs(TRUE)
R <- if (length(args)) as.integer(args[1]) else 200L
.libPaths(c("C:/Users/adf44/source/r/wt-postfit2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({
  library(brms)
  library(frmtmb)
  library(frmtmb.sample)
})
options(mc.cores = 1)
say <- function(...) cat(sprintf(...), "\n", sep = "")
set.seed(9)
dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
fit <- frm(frmtmb::bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
nd <- 100
tpl <- fit$frame$par_template
est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
set.seed(1)
M <- matrix(rep(est, each = nd) + rnorm(nd * length(est), 0, 0.05), nd,
            dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
N <- frmtmb.sample:::draws_to_natural(M, fit)
ds <- structure(list(stanfit = NULL, draws = cbind(N, lp__ = 0), fit = fit),
                class = "frmtmb_draws")
inits <- lapply(seq_len(nd), function(i) {
  s1 <- exp(N[i, "theta_1"])
  list(b = array(N[i, "b_x"], 1),
       Intercept = N[i, "b_Intercept"] + mean(dd$x) * N[i, "b_x"],
       sd_1 = array(s1, 1),
       z_1 = matrix(N[i, sprintf("r_g[%d,Intercept]", 1:6)] / s1, 1),
       sigma = N[i, "sigma"])
})
bk <- suppressMessages(suppressWarnings(
  brm(y ~ x + (1 | g), data = dd, algorithm = "fixed_param", chains = nd,
      iter = 1, warmup = 0, init = inits, refresh = 0, seed = 1,
      silent = 2)))
sd_b <- as.numeric(as_draws_matrix(bk, variable = "sd_g__Intercept"))
sd_f <- exp(N[, "theta_1"])
say("sd per draw, brms vs frmtmb: max rel diff %.3g",
    max(abs(sd_b - sd_f) / sd_f))
spag <- function(o, ...) {
  s <- attr(suppressMessages(suppressWarnings(conditional_effects(
    o, "x", resolution = 3, spaghetti = TRUE, ...)))[[1]], "spaghetti")
  s$estimate__[seq(2, nrow(s), by = 3)]   # the middle grid point
}
pop_b <- spag(bk)
pop_f <- spag(ds)
say("population spaghetti, brms vs frmtmb: max diff %.3g",
    max(abs(pop_b - pop_f)))
zb <- zf <- numeric(0)
wb <- wf <- numeric(R)
for (r in seq_len(R)) {
  set.seed(r)
  ub <- spag(bk, re_formula = NULL, sample_new_levels = "gaussian") - pop_b
  uf <- spag(ds, re_formula = NULL, seed = r) - pop_f
  zb <- c(zb, ub / sd_b)
  zf <- c(zf, uf / sd_f)
  # the band width at the middle point, as each package reports it
  wb[r] <- diff(stats::quantile(pop_b + ub, c(0.025, 0.975)))
  wf[r] <- diff(stats::quantile(pop_f + uf, c(0.025, 0.975)))
}
se_var <- sqrt(2 / length(zb))
say("max |z_brms - z_frmtmb| at matched seeds: %.3g", max(abs(zb - zf)))
say("R = %d repeats x %d draws = %d new-level effects per package", R, nd,
    length(zb))
say("var(z): brms %.4f, frmtmb %.4f (SE of each about %.4f)", var(zb),
    var(zf), se_var)
say("mean(z): brms %.4f, frmtmb %.4f (SE about %.4f)", mean(zb), mean(zf),
    1 / sqrt(length(zb)))
say("quantiles of z (1%%, 10%%, 50%%, 90%%, 99%%): brms %s | frmtmb %s",
    paste(sprintf("%.3f", quantile(zb, c(.01, .1, .5, .9, .99))), collapse = " "),
    paste(sprintf("%.3f", quantile(zf, c(.01, .1, .5, .9, .99))), collapse = " "))
say("ks.test brms z vs frmtmb z: p = %.3f", ks.test(zb, zf)$p.value)
say("middle-point band width over %d repeats: brms %.4f (sd %.4f), frmtmb %.4f (sd %.4f); ratio %.4f, SE of the ratio about %.4f",
    R, mean(wb), sd(wb), mean(wf), sd(wf), mean(wf) / mean(wb),
    (mean(wf) / mean(wb)) * sqrt(var(wb) / mean(wb)^2 / R + var(wf) / mean(wf)^2 / R))
# the reported band itself, frmtmb's own summary, for 30 repeats
bw <- fw <- numeric(min(30L, R))
for (r in seq_len(min(30L, R))) {
  set.seed(1000 + r)
  a <- suppressWarnings(conditional_effects(bk, "x", resolution = 3,
                                            re_formula = NULL,
                                            sample_new_levels = "gaussian"))[[1]]
  b <- suppressMessages(conditional_effects(ds, "x", resolution = 3,
                                            re_formula = NULL,
                                            seed = 1000 + r))[[1]]
  bw[r] <- mean(a$upper__ - a$lower__)
  fw[r] <- mean(b$upper__ - b$lower__)
}
say("reported band width, 30 repeats: brms %.4f (sd %.4f), frmtmb %.4f (sd %.4f), ratio %.4f",
    mean(bw), sd(bw), mean(fw), sd(fw), mean(fw) / mean(bw))
say("done")
