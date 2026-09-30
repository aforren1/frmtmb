# Reviewer: conditional_effects(re_formula = NULL, conditions = ...) on
# draws against brms 2.23.0 at the same draws, on shapes the lane did
# not compare: crossed groups (both set; one set with the other's
# effects pinned at zero so brms's random new-level choice is exact),
# a distributional sigma ~ (1 | g), an ordinal categorical display, an
# unseen level and rows mixing a level with NA. brms is sampled with
# algorithm = "fixed_param", one chain per draw, initialized at the
# draw (the lane's construction in dev/postfit2-celevel-brms.R).
# Draws: the ML estimate plus N(0, 0.05) noise, seed 1.
#   Rscript dev/postfit2-rev-brms.R <base|lane>
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "lane"
libs <- c("C:/Users/adf44/source/r/wt-postfit2-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({
  library(brms)
  library(frmtmb)
  library(frmtmb.sample)
})
options(mc.cores = 1)
say <- function(...) cat(sprintf(...), "\n", sep = "")
say("ARM %s: frmtmb.sample from %s", arm, find.package("frmtmb.sample"))
f3 <- function(v) paste(sprintf("%.4f", v), collapse = " ")
arr <- function(v) array(v, dim = length(v))
mx <- function(a, b) max(abs(as.numeric(a) - as.numeric(b)))
hand <- function(fit, n, seed = 1, edit = NULL) {
  tpl <- fit$frame$par_template
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  set.seed(seed)
  M <- matrix(rep(est, each = n) + rnorm(n * length(est), 0, 0.05), n,
              dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
  if (!is.null(edit)) M <- edit(M)
  list(raw = M, ds = structure(
    list(stanfit = NULL,
         draws = cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0),
         fit = fit), class = "frmtmb_draws"))
}
bfix <- function(formula, data, inits, family = gaussian()) {
  suppressMessages(suppressWarnings(
    brm(formula, data = data, family = family, algorithm = "fixed_param",
        chains = length(inits), iter = 1, warmup = 0, init = inits,
        refresh = 0, seed = 1, silent = 2)))
}
cmp <- function(label, bk, ds, ...) {
  a <- tryCatch(suppressWarnings(conditional_effects(bk, ...)),
                error = function(e) e)
  b <- tryCatch(suppressMessages(conditional_effects(ds, ...)),
                error = function(e) e)
  if (inherits(a, "error") || inherits(b, "error")) {
    say("%s: brms %s | frmtmb %s", label,
        if (inherits(a, "error")) paste("ERROR", conditionMessage(a)) else "ok",
        if (inherits(b, "error")) paste("ERROR", conditionMessage(b)) else "ok")
    return(invisible(NULL))
  }
  a <- a[[1]]
  b <- b[[1]]
  say("%s: rows %d/%d; max |brms - frmtmb| estimate %.3g se %.3g lower %.3g upper %.3g (scale %.3g)",
      label, nrow(a), nrow(b), mx(a$estimate__, b$estimate__),
      mx(a$se__, b$se__), mx(a$lower__, b$lower__),
      mx(a$upper__, b$upper__), max(abs(a$estimate__)))
  invisible(list(a = a, b = b))
}

## 1. crossed (1 | g) + (1 | h), data seed 22
set.seed(22)
db <- data.frame(x = rnorm(160), g = factor(rep(1:8, 20)),
                 h = factor(rep(1:10, each = 16)))
db$y <- rnorm(160, 1 + 0.5 * db$x + rnorm(8, 0, 2)[db$g] +
                rnorm(10, 0, 0.5)[db$h], 1)
fb <- frm(bf(y ~ x + (1 | g) + (1 | h)), family = gaussian(), data = db)
crossed_inits <- function(M) {
  lapply(seq_len(nrow(M)), function(i) {
    s1 <- exp(M[i, "theta_1"])
    s2 <- exp(M[i, "theta_2"])
    list(b = arr(M[i, "b_x"]),
         Intercept = M[i, "b_Intercept"] + mean(db$x) * M[i, "b_x"],
         sd_1 = arr(s1), z_1 = matrix(M[i, sprintf("r_g[%d,Intercept]", 1:8)] / s1, 1),
         sd_2 = arr(s2), z_2 = matrix(M[i, sprintf("r_h[%d,Intercept]", 1:10)] / s2, 1),
         sigma = M[i, "sigma"])
  })
}
h1 <- hand(fb, 5)
bk1 <- bfix(y ~ x + (1 | g) + (1 | h), db, crossed_inits(h1$ds$draws))
say("1 crossed: brms r_g / r_h vs frmtmb: %.3g %.3g",
    mx(as_draws_matrix(bk1, variable = "r_g[3,Intercept]"),
       h1$ds$draws[, "r_g[3,Intercept]"]),
    mx(as_draws_matrix(bk1, variable = "r_h[4,Intercept]"),
       h1$ds$draws[, "r_h[4,Intercept]"]))
cmp("1a both set g=7,h=4", bk1, h1$ds, effects = "x", re_formula = NULL,
    resolution = 5,
    conditions = data.frame(g = factor(7, levels = 1:8),
                            h = factor(4, levels = 1:10)))
## h's effects pinned at 0 and its sd at 1e-8: brms's new level picks an
## existing level's draw, 0; frmtmb draws N(0, 1e-8)
h2 <- hand(fb, 5, edit = function(M) {
  M[, grep("^r_h", colnames(M))] <- 0
  M[, "theta_2"] <- log(1e-8)
  M
})
bk2 <- bfix(y ~ x + (1 | g) + (1 | h), db, crossed_inits(h2$ds$draws))
cmp("1b g=7 set, h unset (h pinned at 0)", bk2, h2$ds, effects = "x",
    re_formula = NULL, resolution = 5, seed = 1,
    conditions = data.frame(g = factor(7, levels = 1:8)))
cmp("1c x:g effect, h unset (h pinned)", bk2, h2$ds, effects = "x:g",
    re_formula = NULL, resolution = 5, seed = 1)
cmp("1d re_formula ~(1|g), g=7", bk1, h1$ds, effects = "x",
    re_formula = ~ (1 | g), resolution = 5,
    conditions = data.frame(g = factor(7, levels = 1:8)))

## 2. distributional sigma ~ (1 | g)
fs <- frm(bf(y ~ x + (1 | g), sigma ~ (1 | g)), family = gaussian(),
          data = db)
hs <- hand(fs, 5)
sig_inits <- lapply(seq_len(5), function(i) {
  M <- hs$ds$draws
  s1 <- exp(M[i, "theta_1"])
  s2 <- exp(M[i, "theta_2"])
  list(b = arr(M[i, "b_x"]),
       Intercept = M[i, "b_Intercept"] + mean(db$x) * M[i, "b_x"],
       Intercept_sigma = M[i, "b_sigma_Intercept"],
       sd_1 = arr(s1), z_1 = matrix(M[i, sprintf("r_g[%d,Intercept]", 1:8)] / s1, 1),
       sd_2 = arr(s2),
       z_2 = matrix(M[i, sprintf("r_g__sigma[%d,Intercept]", 1:8)] / s2, 1))
})
bks <- bfix(brms::bf(y ~ x + (1 | g), sigma ~ (1 | g)), db, sig_inits)
sv <- grep("^r_g", variables(bks), value = TRUE)
say("2 sigma: brms group variables %s", paste(head(sv, 3), tail(sv, 3),
                                              collapse = " "))
bsig <- grep("sigma", sv, value = TRUE)[3]
say("2 sigma: brms %s vs frmtmb r_g__sigma[3,Intercept] %.3g", bsig,
    mx(as_draws_matrix(bks, variable = bsig),
       hs$ds$draws[, "r_g__sigma[3,Intercept]"]))
cmp("2a mu, g=3", bks, hs$ds, effects = "x", re_formula = NULL,
    resolution = 5, conditions = data.frame(g = factor(3, levels = 1:8)))
cmp("2b dpar sigma, g=3", bks, hs$ds, effects = "x", re_formula = NULL,
    resolution = 5, dpar = "sigma",
    conditions = data.frame(g = factor(3, levels = 1:8)))

## 3. unseen level and mixed NA rows, (1 | g); 40 draws. brms draws a
## new level from a random existing level per draw, so only the widths
## are compared
set.seed(9)
dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
f1 <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
h3 <- hand(f1, 40)
g_inits <- lapply(seq_len(40), function(i) {
  M <- h3$ds$draws
  s1 <- exp(M[i, "theta_1"])
  list(b = arr(M[i, "b_x"]),
       Intercept = M[i, "b_Intercept"] + mean(dd$x) * M[i, "b_x"],
       sd_1 = arr(s1),
       z_1 = matrix(M[i, sprintf("r_g[%d,Intercept]", 1:6)] / s1, 1),
       sigma = M[i, "sigma"])
})
bk3 <- bfix(y ~ x + (1 | g), dd, g_inits)
width <- function(d) f3(d$upper__ - d$lower__)
wb <- function(o, ...) tryCatch(suppressMessages(suppressWarnings(
  conditional_effects(o, effects = "x", resolution = 3, seed = 1, ...))),
  error = function(e) e)
for (case in list(
  list(nm = "pop (re_formula NA)", a = list()),
  list(nm = "observed g=2", a = list(re_formula = NULL,
    conditions = data.frame(g = factor(2, levels = 1:6)))),
  list(nm = "unseen g='99'", a = list(re_formula = NULL,
    conditions = data.frame(g = "99"))),
  list(nm = "unseen g='99' allow_new_levels", a = list(re_formula = NULL,
    conditions = data.frame(g = "99"), allow_new_levels = TRUE)),
  list(nm = "unset (no conditions)", a = list(re_formula = NULL)))) {
  a <- do.call(wb, c(list(bk3), case$a))
  b <- do.call(wb, c(list(h3$ds), case$a))
  say("3 %s: width brms %s | frmtmb %s", case$nm,
      if (inherits(a, "error")) paste("ERROR", substr(conditionMessage(a), 1, 80)) else width(a[[1]]),
      if (inherits(b, "error")) paste("ERROR", substr(conditionMessage(b), 1, 80)) else width(b[[1]]))
}
cm <- data.frame(g = factor(c("2", NA), levels = 1:6))
rownames(cm) <- c("lev2", "unset")
for (anl in c(FALSE, TRUE)) {
  a <- wb(bk3, re_formula = NULL, conditions = cm)  # brms sets allow_new_levels = TRUE itself
  b <- wb(h3$ds, re_formula = NULL, conditions = cm, allow_new_levels = anl)
  for (cv in c("lev2", "unset")) {
    say("3 mixed rows anl=%s cond %s: width brms %s | frmtmb %s", anl, cv,
        if (inherits(a, "error")) paste("ERROR", substr(conditionMessage(a), 1, 80)) else
          width(a[[1]][a[[1]]$cond__ == cv, ]),
        if (inherits(b, "error")) paste("ERROR", substr(conditionMessage(b), 1, 80)) else
          width(b[[1]][b[[1]]$cond__ == cv, ]))
  }
}
say("done")
