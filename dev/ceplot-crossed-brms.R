# Lane ceplot: crossed (1 | g) + (1 | h) + (1 | g:h) at g and h set to
# observed levels never seen together, on draws, against brms 2.23.0 at
# the same draws with sample_new_levels = "gaussian"; and on the fit,
# the bootstrap band against the Wald band.
# brms is sampled with algorithm = "fixed_param", one chain per draw,
# initialized at the draw (the construction of dev/postfit2-p2-brms.R).
# Both packages are called with the same seed (set.seed(s) before brms,
# seed = s for frmtmb): drawing the same random numbers in the same
# order, the bands agree to rounding.
#   Rscript dev/ceplot-crossed-brms.R > dev/ceplot-log/crossed-brms.txt
# Seeds: data 49 (dev/postfit2-p2-checks.R's P1-M2), draws 1, curves 1
# to 3, bootstrap 5.
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(brms)
  library(frmtmb)
  library(frmtmb.sample)
})
options(mc.cores = 1)
cat("frmtmb from", find.package("frmtmb"), "\n")
say <- function(...) cat(sprintf(...), "\n", sep = "")
f4 <- function(v) paste(sprintf("%.4f", v), collapse = " ")
mx <- function(a, b) max(abs(a - b))
arr <- function(v) array(v, dim = length(v))
N <- 200
hand <- function(fit, n, seed = 1) {
  tpl <- fit$frame$par_template
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  set.seed(seed)
  M <- matrix(rep(est, each = n) + rnorm(n * length(est), 0, 0.05), n,
              dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
  structure(list(stanfit = NULL,
                 draws = cbind(frmtmb.sample:::draws_to_natural(M, fit),
                               lp__ = 0),
                 fit = fit), class = "frmtmb_draws")
}
set.seed(49)
dc <- expand.grid(g = factor(1:6), h = factor(1:5), r = 1:4)
dc <- dc[!(dc$g == "1" & dc$h == "1"), ]
dc$x <- rnorm(nrow(dc))
dc$y <- rnorm(nrow(dc), 1 + 0.5 * dc$x + rnorm(6)[dc$g] + rnorm(5)[dc$h] +
                rnorm(30, 0, 0.7)[as.integer(interaction(dc$g, dc$h))], 0.5)
fc <- frm(bf(y ~ x + (1 | g) + (1 | h) + (1 | g:h)), family = gaussian(),
          data = dc)
vc <- VarCorr(fc)
say("fit sd: g %.4f h %.4f g:h %.4f", vc$g$sd[1, 1], vc$h$sd[1, 1],
    vc$`g:h`$sd[1, 1])
ds <- hand(fc, N)
M <- ds$draws
# theta_k is the log sd of the k-th block, blocks in the labels' order
th <- exp(M[, c("theta_1", "theta_2", "theta_3")])
say("theta order check (block sd at the estimate): %s",
    f4(exp(fc$estimates$theta)))

## brms at the same draws: sd_1 = g, sd_2 = g:h, sd_3 = h
## (dev/ceplot-log/crossed-map.txt)
lab_g <- sprintf("r_g[%d,Intercept]", 1:6)
lab_h <- sprintf("r_h[%d,Intercept]", 1:5)
lab_gh <- grep("^r_g:h\\[", colnames(M), value = TRUE)
inits <- lapply(seq_len(N), function(i) {
  list(b = arr(M[i, "b_x"]),
       Intercept = M[i, "b_Intercept"] + mean(dc$x) * M[i, "b_x"],
       sigma = M[i, "sigma"],
       sd_1 = arr(th[i, 1]), sd_2 = arr(th[i, 3]), sd_3 = arr(th[i, 2]),
       z_1 = matrix(M[i, lab_g] / th[i, 1], 1),
       z_2 = matrix(M[i, lab_gh] / th[i, 3], 1),
       z_3 = matrix(M[i, lab_h] / th[i, 2], 1))
})
b <- suppressMessages(suppressWarnings(
  brm(y ~ x + (1 | g) + (1 | h) + (1 | g:h), data = dc,
      algorithm = "fixed_param", chains = N, iter = 1, warmup = 0,
      init = inits, refresh = 0, seed = 1, silent = 2)))
bm <- as_draws_matrix(b)
say("brms draws %d; r_g:h max diff %.3g; sd_g:h max rel diff %.3g",
    ndraws(b), mx(bm[, lab_gh], M[, lab_gh]),
    max(abs(bm[, "sd_g:h__Intercept"] / th[, 3] - 1)))
ce <- function(o, cond, ...) {
  suppressMessages(suppressWarnings(conditional_effects(
    o, effects = "x", resolution = 3, re_formula = NULL,
    conditions = cond, ...)))[[1]]
}
obs <- data.frame(g = "1", h = "2")
a <- ce(b, obs)
d <- ce(ds, obs)
say("observed g = 1, h = 2 (no draw): max diff estimate %.3g lower %.3g upper %.3g",
    mx(a$estimate__, d$estimate__), mx(a$lower__, d$lower__),
    mx(a$upper__, d$upper__))
new <- data.frame(g = "1", h = "1")
for (s in 1:3) {
  set.seed(s)
  a <- ce(b, new, sample_new_levels = "gaussian")
  d <- tryCatch(ce(ds, new, seed = s), error = function(e) e)
  if (inherits(d, "error")) {
    say("unseen g:h, seed %d: frmtmb ERROR %s", s, conditionMessage(d))
    next
  }
  say("unseen g:h, seed %d: estimate brms %s | frmtmb %s", s,
      f4(a$estimate__), f4(d$estimate__))
  say("unseen g:h, seed %d: lower brms %s | frmtmb %s; upper brms %s | frmtmb %s",
      s, f4(a$lower__), f4(d$lower__), f4(a$upper__), f4(d$upper__))
  say("unseen g:h, seed %d: max diff estimate %.3g lower %.3g upper %.3g", s,
      mx(a$estimate__, d$estimate__), mx(a$lower__, d$lower__),
      mx(a$upper__, d$upper__))
}

## the fit: bootstrap band against the Wald band
w_new <- ce(fc, list(g = "1", h = "1"))
w_obs <- ce(fc, list(g = "1", h = "2"))
bo <- tryCatch(ce(fc, list(g = "1", h = "1"), band = "boot", boot = 60,
                  seed = 5), error = function(e) e)
say("fit unseen g:h | wald: est %s | width %s", f4(w_new$estimate__),
    f4(w_new$upper__ - w_new$lower__))
say("fit observed g:h | wald: est %s | width %s", f4(w_obs$estimate__),
    f4(w_obs$upper__ - w_obs$lower__))
if (inherits(bo, "error")) {
  say("fit unseen g:h | boot60: ERROR %s", conditionMessage(bo))
} else {
  say("fit unseen g:h | boot60: est %s | width %s", f4(bo$estimate__),
      f4(bo$upper__ - bo$lower__))
  say("fit unseen g:h: boot width / wald width %s",
      f4((bo$upper__ - bo$lower__) / (w_new$upper__ - w_new$lower__)))
}
dr <- ce(ds, new, seed = 1)
do <- ce(ds, obs, seed = 1)
w_new_w <- w_new$upper__ - w_new$lower__
say("draws unseen g:h: width %s (%s of the Wald width); draws observed g:h: width %s",
    f4(dr$upper__ - dr$lower__), f4((dr$upper__ - dr$lower__) / w_new_w),
    f4(do$upper__ - do$lower__))
say("done")
