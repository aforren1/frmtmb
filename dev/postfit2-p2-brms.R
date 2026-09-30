# Punch round 2: gr(g, by = f) and mm(g1, g2) on draws against brms
# 2.23.0 at the same draws, brms with sample_new_levels = "gaussian".
# brms is sampled with algorithm = "fixed_param", one chain per draw,
# initialized at the draw (as dev/postfit2-p1-brms.R does). Where a new
# level is drawn, both packages are called with the same seed
# (set.seed(s) before brms, seed = s for frmtmb): if they draw the same
# random numbers in the same order, the bands agree to rounding.
#   Rscript dev/postfit2-p2-brms.R > dev/postfit2-log/p2-brms.txt
# Seeds: data 44 (gr by) and 43 (mm), draws 1, curves 1 to 3.
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
bfix <- function(formula, data, inits) {
  suppressMessages(suppressWarnings(
    brm(formula, data = data, algorithm = "fixed_param",
        chains = length(inits), iter = 1, warmup = 0, init = inits,
        refresh = 0, seed = 1, silent = 2)))
}
ce <- function(o, ...) {
  suppressMessages(suppressWarnings(conditional_effects(
    o, effects = "x", resolution = 3, re_formula = NULL, ...)))[[1]]
}
cmp <- function(label, cond, bo, fo, seeds = 1:3) {
  for (s in seeds) {
    set.seed(s)
    a <- tryCatch(ce(bo, conditions = cond, sample_new_levels = "gaussian"),
                  error = function(e) e)
    b <- tryCatch(ce(fo, conditions = if (is.null(cond)) list() else cond,
                     seed = s), error = function(e) e)
    if (inherits(a, "error") || inherits(b, "error")) {
      say("%s, seed %d: brms %s | frmtmb %s", label, s,
          if (inherits(a, "error")) paste("ERROR", conditionMessage(a))
          else "answers",
          if (inherits(b, "error")) paste("ERROR", conditionMessage(b))
          else "answers")
      next
    }
    say("%s, seed %d: estimate brms %s | frmtmb %s", label, s,
        f4(a$estimate__), f4(b$estimate__))
    say("%s, seed %d: lower brms %s | frmtmb %s; upper brms %s | frmtmb %s",
        label, s, f4(a$lower__), f4(b$lower__), f4(a$upper__),
        f4(b$upper__))
    say("%s, seed %d: max diff estimate %.3g lower %.3g upper %.3g", label,
        s, mx(a$estimate__, b$estimate__), mx(a$lower__, b$lower__),
        mx(a$upper__, b$upper__))
  }
}

## 1. (1 | gr(g, by = f)), data seed 44
set.seed(44)
d4 <- data.frame(x = rnorm(240), g = factor(rep(1:12, 20)))
d4$f <- factor(ifelse(as.integer(d4$g) <= 6, "a", "b"))
sdg <- ifelse(1:12 <= 6, 0.3, 1.5)
d4$y <- rnorm(240, 1 + 0.5 * d4$x + rnorm(12, 0, sdg)[d4$g], 0.5)
fit4 <- frm(bf(y ~ x + f + (1 | gr(g, by = f))), family = gaussian(),
            data = d4)
ds4 <- hand(fit4, N)
M <- ds4$draws
mfb <- mean(d4$f == "b")
by_of <- as.character(d4$f[match(1:12, as.integer(d4$g))])
inits4 <- lapply(seq_len(N), function(i) {
  sds <- exp(M[i, c("theta_1", "theta_2")])
  r <- M[i, sprintf("r_g[%d,Intercept]", 1:12)]
  list(b = arr(M[i, c("b_x", "b_fb")]),
       Intercept = M[i, "b_Intercept"] + mean(d4$x) * M[i, "b_x"] +
         mfb * M[i, "b_fb"],
       sigma = M[i, "sigma"],
       sd_1 = matrix(sds, 1),
       z_1 = matrix(r / sds[ifelse(by_of == "a", 1, 2)], 1))
})
b4 <- bfix(y ~ x + f + (1 | gr(g, by = f)), d4, inits4)
say("1: draws brms %d; r_g[3] max diff %.3g; r_g[9] max diff %.3g",
    ndraws(b4),
    mx(as_draws_matrix(b4, variable = "r_g[3,Intercept]"),
       M[, "r_g[3,Intercept]"]),
    mx(as_draws_matrix(b4, variable = "r_g[9,Intercept]"),
       M[, "r_g[9,Intercept]"]))
say("1: sd_g__Intercept:fa max rel diff %.3g",
    max(abs(as_draws_matrix(b4, variable = "sd_g__Intercept:fa") /
              exp(M[, "theta_1"]) - 1)))
for (cond in list(data.frame(f = "a", g = "3"), data.frame(f = "b", g = "9"))) {
  a <- ce(b4, conditions = cond)
  b <- ce(ds4, conditions = cond)
  say("1 observed f = %s, g = %s (no draw): max diff estimate %.3g lower %.3g upper %.3g",
      cond$f, cond$g, mx(a$estimate__, b$estimate__),
      mx(a$lower__, b$lower__), mx(a$upper__, b$upper__))
}
cmp("1 f = a, g unset", data.frame(f = "a"), b4, ds4)
cmp("1 f = b, g unset", data.frame(f = "b"), b4, ds4)
cmp("1 f = a, g = 99", data.frame(f = "a", g = "99"), b4, ds4)
cmp("1 f = b, g = 99", data.frame(f = "b", g = "99"), b4, ds4)

## 2. (1 | mm(g1, g2)), data seed 43
set.seed(43)
d3 <- data.frame(x = rnorm(200), g1 = factor(sample(1:10, 200, TRUE)),
                 g2 = factor(sample(1:10, 200, TRUE)))
u <- rnorm(10, 0, 1)
d3$y <- rnorm(200, 1 + 0.5 * d3$x + 0.5 * (u[d3$g1] + u[d3$g2]), 0.5)
fit3 <- frm(bf(y ~ x + (1 | mm(g1, g2))), family = gaussian(), data = d3)
ds3 <- hand(fit3, N)
M <- ds3$draws
inits3 <- lapply(seq_len(N), function(i) {
  s1 <- exp(M[i, "theta_1"])
  list(b = arr(M[i, "b_x"]),
       Intercept = M[i, "b_Intercept"] + mean(d3$x) * M[i, "b_x"],
       sigma = M[i, "sigma"], sd_1 = arr(s1),
       z_1 = matrix(M[i, sprintf("r_mmg1g2[%d,Intercept]", 1:10)] / s1, 1))
})
b3 <- bfix(y ~ x + (1 | mm(g1, g2)), d3, inits3)
say("2: draws brms %d; r_mmg1g2[2] max diff %.3g", ndraws(b3),
    mx(as_draws_matrix(b3, variable = "r_mmg1g2[2,Intercept]"),
       M[, "r_mmg1g2[2,Intercept]"]))
cond <- data.frame(g1 = "2", g2 = "2")
a <- ce(b3, conditions = cond)
b <- ce(ds3, conditions = cond)
say("2 g1 = g2 = 2 (no draw): max diff estimate %.3g lower %.3g upper %.3g",
    mx(a$estimate__, b$estimate__), mx(a$lower__, b$lower__),
    mx(a$upper__, b$upper__))
cmp("2 nothing set", NULL, b3, ds3)
cmp("2 g1 = 2, g2 unset", data.frame(g1 = "2"), b3, ds3)
cmp("2 g1 = 99, g2 = 98", data.frame(g1 = "99", g2 = "98"), b3, ds3)
say("done")
