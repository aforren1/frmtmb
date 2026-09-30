# Punch round 1: new levels on draws against brms 2.23.0 at the same
# draws, with brms's sample_new_levels = "gaussian", the law frmtmb
# draws a new level from. A new level is a random draw in both, so the
# comparison there is of band WIDTHS over 400 draws; where no new level
# is drawn the values are compared exactly. brms is sampled with
# algorithm = "fixed_param", one chain per draw, initialized at the
# draw; the r_ columns are matched to brms's by name, which frmtmb's
# draws labels mirror, and each group's sd is exp() of its theta.
#   Rscript dev/postfit2-p1-brms.R > dev/postfit2-log/p1-brms.txt
# Seeds: data 9 (one term) and 25 (nested), draws 1, curves 1.
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
f4 <- function(v) paste(sprintf("%.4f", v), collapse = " ")
arr <- function(v) array(v, dim = length(v))
mx <- function(a, b) max(abs(a - b))
N <- 400
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
# `theta`: brms group name -> frmtmb theta column, in brms's group order
inits_for <- function(ds, bv, data, theta) {
  M <- ds$draws
  lapply(seq_len(nrow(M)), function(i) {
    out <- list(b = arr(M[i, "b_x"]),
                Intercept = M[i, "b_Intercept"] + mean(data$x) * M[i, "b_x"],
                sigma = M[i, "sigma"])
    for (k in seq_along(theta)) {
      rn <- bv[startsWith(bv, paste0("r_", names(theta)[k], "["))]
      sdk <- exp(M[i, theta[[k]]])
      out[[paste0("sd_", k)]] <- arr(sdk)
      out[[paste0("z_", k)]] <- matrix(M[i, rn] / sdk, 1)
    }
    out
  })
}
bfix <- function(formula, data, inits, fit = NA) {
  suppressMessages(suppressWarnings(
    brm(formula, data = data, algorithm = "fixed_param", fit = fit,
        chains = length(inits), iter = 1, warmup = 0, init = inits,
        refresh = 0, seed = 1, silent = 2)))
}
ce <- function(o, ...) {
  suppressMessages(suppressWarnings(conditional_effects(
    o, effects = "x", resolution = 3, re_formula = NULL, ...)))[[1]]
}
wd <- function(d) d$upper__ - d$lower__
report <- function(label, a, b) {
  say("%s: width brms %s | frmtmb %s | ratio %s", label, f4(wd(a)),
      f4(wd(b)), f4(wd(b) / wd(a)))
}

## 1. (1 | g), data seed 9
set.seed(9)
dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
f1 <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
d1 <- hand(f1, N)
b0 <- suppressMessages(suppressWarnings(
  brm(y ~ x + (1 | g), data = dd, algorithm = "fixed_param", chains = 1,
      iter = 1, warmup = 0, refresh = 0, seed = 1, silent = 2)))
bv <- variables(b0)
b1 <- bfix(y ~ x + (1 | g), dd, inits_for(d1, bv, dd, c(g = "theta_1")),
           fit = b0)
say("1: draws brms %d; r_g[2] max diff %.3g", ndraws(b1),
    mx(as_draws_matrix(b1, variable = "r_g[2,Intercept]"),
       d1$draws[, "r_g[2,Intercept]"]))
g2 <- list(g = factor("2", levels = 1:6))
a <- ce(b1, conditions = g2)
b <- ce(d1, conditions = g2)
say("1 observed g=2 (no draw): estimate %.3g lower %.3g upper %.3g",
    mx(a$estimate__, b$estimate__), mx(a$lower__, b$lower__),
    mx(a$upper__, b$upper__))
pop <- suppressMessages(conditional_effects(d1, effects = "x", resolution = 3))[[1]]
report("1 population, for scale", pop, pop)
for (cs in list(list(nm = "unset g (no conditions)", c = list()),
                list(nm = "unseen g = '99'", c = list(g = "99")))) {
  report(paste("1", cs$nm), ce(b1, conditions = if (length(cs$c)) cs$c,
                               sample_new_levels = "gaussian"),
         ce(d1, conditions = cs$c, seed = 1))
}
cm <- data.frame(g = factor(c("2", NA), levels = 1:6))
rownames(cm) <- c("lev2", "unset")
a <- ce(b1, conditions = cm, sample_new_levels = "gaussian")
b <- ce(d1, conditions = cm, seed = 1)
say("1 mixed rows, lev2 (no draw): upper max diff %.3g",
    mx(a$upper__[a$cond__ == "lev2"], b$upper__[b$cond__ == "lev2"]))
report("1 mixed rows, unset", a[a$cond__ == "unset", ],
       b[b$cond__ == "unset", ])

# the widths above are one random realization each: averaged over 8
# brms calls and 8 frmtmb seeds, the Monte Carlo spread shows beside it
rep_w <- function(f) do.call(rbind, lapply(1:8, f))
wb8 <- rep_w(function(s) wd(ce(b1, sample_new_levels = "gaussian")))
wf8 <- rep_w(function(s) wd(ce(d1, conditions = list(), seed = s)))
say("1 unset g, mean width over 8: brms %s (sd %s) | frmtmb %s (sd %s) | ratio %s",
    f4(colMeans(wb8)), f4(apply(wb8, 2, sd)), f4(colMeans(wf8)),
    f4(apply(wf8, 2, sd)), f4(colMeans(wf8) / colMeans(wb8)))

## 2. nested (1 | g / h), data seed 25
set.seed(25)
dn <- expand.grid(g = factor(1:8), h = factor(1:4), r = 1:5)
dn$x <- rnorm(nrow(dn))
gh <- as.integer(interaction(dn$g, dn$h))
dn$y <- rnorm(nrow(dn), 1 + 0.5 * dn$x + rnorm(8, 0, 2)[dn$g] +
                rnorm(32, 0, 0.7)[gh])
f2 <- frm(bf(y ~ x + (1 | g / h)), family = gaussian(), data = dn)
d2 <- hand(f2, N)
b0n <- suppressMessages(suppressWarnings(
  brm(y ~ x + (1 | g / h), data = dn, algorithm = "fixed_param",
      chains = 1, iter = 1, warmup = 0, refresh = 0, seed = 1,
      silent = 2)))
bvn <- variables(b0n)
grp <- sub("^sd_(.*)__Intercept$", "\\1",
           grep("^sd_.*__Intercept$", bvn, value = TRUE))
say("2: brms groups in order: %s", paste(grp, collapse = ", "))
# frmtmb's blocks: "h:g" holds theta_1, "g" theta_2
theta <- c(g = "theta_2", "g:h" = "theta_1")[grp]
b2 <- bfix(y ~ x + (1 | g / h), dn, inits_for(d2, bvn, dn, theta),
           fit = b0n)
say("2: r_g:h[5_2] max diff %.3g",
    mx(as_draws_matrix(b2, variable = "r_g:h[5_2,Intercept]"),
       d2$draws[, "r_g:h[5_2,Intercept]"]))
both <- list(g = factor("5", levels = 1:8), h = factor("2", levels = 1:4))
a <- ce(b2, conditions = both)
b <- ce(d2, conditions = both)
say("2 g=5, h=2 both observed (no draw): estimate %.3g upper %.3g",
    mx(a$estimate__, b$estimate__), mx(a$upper__, b$upper__))
g5 <- list(g = factor("5", levels = 1:8))
a <- ce(b2, conditions = g5, sample_new_levels = "gaussian")
b <- ce(d2, conditions = g5, seed = 1)
report("2 g=5 set, h unset (new g:h within g = 5)", a, b)
say("2 g=5 set, h unset: estimate brms %s | frmtmb %s", f4(a$estimate__),
    f4(b$estimate__))
say("done")
