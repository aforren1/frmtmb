# Copy of dev/postfit2-rev-p2-mm.R for the ceplot review, with the library
# paths of this round (wt-ceplot-lib, rellib-r4); otherwise unchanged.
# Reviewer, punch round 2: mm() on draws against brms 2.23.0 at the same
# draws (brms fixed_param, one chain per draw; sample_new_levels =
# "gaussian"; matched seeds s = 1:2 as in dev/postfit2-p2-brms.R), on
# shapes the lane did not compare: unequal weights, three members, a
# repeated unseen level, mm() beside a plain (1 | g1); and the Wald
# band's new-level variance at a mixed row, against w^2 * sd^2.
#   Rscript dev/postfit2-rev-p2-mm.R
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
say <- function(...) cat(sprintf(...), "\n", sep = "")
f4 <- function(v) paste(sprintf("%.4f", v), collapse = " ")
arr <- function(v) array(v, dim = length(v))
mx <- function(a, b) max(abs(a - b))
N <- 100
hand <- function(fit, n, seed = 1) {
  tpl <- fit$frame$par_template
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  set.seed(seed)
  M <- matrix(rep(est, each = n) + rnorm(n * length(est), 0, 0.05), n,
              dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
  structure(list(stanfit = NULL,
                 draws = cbind(frmtmb.sample:::draws_to_natural(M, fit),
                               lp__ = 0), fit = fit), class = "frmtmb_draws")
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
cmp <- function(label, cond, bo, fo, seeds = 1:2) {
  for (s in seeds) {
    set.seed(s)
    a <- tryCatch(ce(bo, conditions = cond, sample_new_levels = "gaussian"),
                  error = function(e) e)
    b <- tryCatch(ce(fo, conditions = if (is.null(cond)) list() else cond,
                     seed = s), error = function(e) e)
    if (inherits(a, "error") || inherits(b, "error")) {
      say("%s, seed %d: brms %s | frmtmb %s", label, s,
          if (inherits(a, "error")) paste("ERROR", substr(conditionMessage(a), 1, 150)) else "answers",
          if (inherits(b, "error")) paste("ERROR", substr(conditionMessage(b), 1, 150)) else "answers")
      next
    }
    say("%s, seed %d: max diff estimate %.3g lower %.3g upper %.3g | width brms %s frmtmb %s",
        label, s, mx(a$estimate__, b$estimate__), mx(a$lower__, b$lower__),
        mx(a$upper__, b$upper__), f4(a$upper__ - a$lower__),
        f4(b$upper__ - b$lower__))
  }
}
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
lev <- sprintf("r_%s[%d,Intercept]", "%s", 1:10)
zmat <- function(M, i, grp, s) {
  matrix(M[i, sprintf("r_%s[%d,Intercept]", grp, 1:10)] / s, 1)
}

## A. unequal weights
fA <- tryCatch(frm(bf(y ~ x + (1 | mm(g1, g2, weights = cbind(w1, w2)))),
                   family = gaussian(), data = d), error = function(e) e)
if (inherits(fA, "error")) say("A fit: ERROR %s", conditionMessage(fA)) else {
  dA <- hand(fA, N)
  M <- dA$draws
  grp <- sub("^r_([^[]+)\\[.*", "\\1", grep("^r_", colnames(M), value = TRUE)[1])
  say("A frmtmb group label %s", grp)
  inA <- lapply(seq_len(N), function(i) {
    s1 <- exp(M[i, "theta_1"])
    list(b = arr(M[i, "b_x"]), Intercept = M[i, "b_Intercept"] +
           mean(d$x) * M[i, "b_x"], sigma = M[i, "sigma"], sd_1 = arr(s1),
         z_1 = zmat(M, i, grp, s1))
  })
  bA <- bfix(y ~ x + (1 | mm(g1, g2, weights = cbind(w1, w2))), d, inA)
  bgrp <- sub("^r_([^[]+)\\[.*", "\\1", grep("^r_", variables(bA), value = TRUE)[1])
  say("A brms group label %s (frmtmb %s)", bgrp, grp)
  say("A draws brms %d; %s max diff %.3g", ndraws(bA), sprintf("r_%s[2,Intercept]", bgrp),
      mx(as_draws_matrix(bA, variable = sprintf("r_%s[2,Intercept]", bgrp)),
         M[, sprintf("r_%s[2,Intercept]", grp)]))
  cmp("A weights, g1 = g2 = 2 (no draw)", data.frame(g1 = "2", g2 = "2"), bA, dA, 1)
  cmp("A weights, nothing set", NULL, bA, dA)
  cmp("A weights, g1 = 2, g2 unset", data.frame(g1 = "2"), bA, dA)
  cmp("A weights, g1 = g2 = 99", data.frame(g1 = "99", g2 = "99"), bA, dA)
  cmp("A weights, g1 = 99, g2 = 98 (brms shares one draw)",
      data.frame(g1 = "99", g2 = "98"), bA, dA)
  ## the Wald band's new-level variance at a mixed row: the extra
  ## variance must be w2^2 * sd^2 at the grid's weights
  w <- ce(fA, conditions = data.frame(g1 = "2"))
  sd1 <- exp(fA$estimates$theta[1])
  w2 <- mean(d$w2)
  lp <- frmtmb:::find_linpred(fA, "y", "mu")
  nd <- w[c("x", "g1", "g2", "w1", "w2")]
  ed <- frmtmb:::lp_eta_design(fA, lp, nd, TRUE, TRUE)
  ev <- tryCatch(frmtmb:::lp_extra_var(fA, ed, TRUE), error = function(e) e)
  if (inherits(ev, "error")) say("A extra var: ERROR %s", conditionMessage(ev)) else {
    num <- unlist(Filter(is.numeric, if (is.list(ev)) ev else list(ev)))
    say("A Wald at g1 = 2, g2 unset: lp_extra_var numeric parts %s; w2^2 sd^2 = %.5f (w1, w2 held at %.4f %.4f)",
        f4(num), w2^2 * sd1^2, mean(d$w1), w2)
  }
  ## the new member's variance in the Wald band: with the fit's sd set
  ## to 1 and to 2 (the coefficient covariance unchanged), se^2 must
  ## differ by w2^2 * (4 - 1) at a mixed row, (w1 + w2)^2 * 3 with
  ## nothing set, and (w1^2 + w2^2) * 3 at two different unseen members
  fA1 <- fA; fA1$estimates$theta[1] <- log(1)
  fA2 <- fA; fA2$estimates$theta[1] <- log(2)
  w1m <- mean(d$w1)
  for (cc in list(list(nm = "g1 = 2, g2 unset", c = data.frame(g1 = "2"),
                       want = 3 * w2^2),
                  list(nm = "nothing set", c = NULL, want = 3 * (w1m + w2)^2),
                  list(nm = "g1 = 99, g2 = 98", c = data.frame(g1 = "99", g2 = "98"),
                       want = 3 * (w1m^2 + w2^2)),
                  list(nm = "g1 = g2 = 99", c = data.frame(g1 = "99", g2 = "99"),
                       want = 3 * (w1m + w2)^2))) {
    cnd <- if (is.null(cc$c)) list() else cc$c
    s1 <- ce(fA1, conditions = cnd)$se__
    s2 <- ce(fA2, conditions = cnd)$se__
    say("A Wald, %s: se^2 difference %s, expected %.5f", cc$nm,
        f4(s2^2 - s1^2), cc$want)
  }
  wo <- ce(fA, conditions = data.frame(g1 = "2", g2 = "2"))
  say("A Wald widths: mixed row %s | both observed %s | nothing set %s",
      f4(w$upper__ - w$lower__), f4(wo$upper__ - wo$lower__),
      f4(ce(fA)$upper__ - ce(fA)$lower__))
}

## B. three members
fB <- frm(bf(y3 ~ x + (1 | mm(g1, g2, g3))), family = gaussian(), data = d)
dB <- hand(fB, N)
M <- dB$draws
grpB <- sub("^r_([^[]+)\\[.*", "\\1", grep("^r_", colnames(M), value = TRUE)[1])
inB <- lapply(seq_len(N), function(i) {
  s1 <- exp(M[i, "theta_1"])
  list(b = arr(M[i, "b_x"]), Intercept = M[i, "b_Intercept"] +
         mean(d$x) * M[i, "b_x"], sigma = M[i, "sigma"], sd_1 = arr(s1),
       z_1 = zmat(M, i, grpB, s1))
})
bB <- bfix(y3 ~ x + (1 | mm(g1, g2, g3)), d, inB)
say("B draws brms %d; r_%s[2] max diff %.3g", ndraws(bB), grpB,
    mx(as_draws_matrix(bB, variable = sprintf("r_%s[2,Intercept]", grpB)),
       M[, sprintf("r_%s[2,Intercept]", grpB)]))
cmp("B three, nothing set", NULL, bB, dB)
cmp("B three, g1 = 2, others unset", data.frame(g1 = "2"), bB, dB)
cmp("B three, g1 = 2, g2 = g3 = 99", data.frame(g1 = "2", g2 = "99", g3 = "99"), bB, dB)
cmp("B three, g1 = 2, g2 = 99, g3 = 98 (brms shares)",
    data.frame(g1 = "2", g2 = "99", g3 = "98"), bB, dB)

## C. mm() beside a plain (1 | g1)
fC <- frm(bf(yp ~ x + (1 | mm(g1, g2)) + (1 | g1)), family = gaussian(),
          data = d)
dC <- hand(fC, N)
M <- dC$draws
say("C frmtmb labels: %s", paste(unique(sub("\\[.*", "", grep("^r_|^theta", colnames(M), value = TRUE))), collapse = " "))
inC <- lapply(seq_len(N), function(i) {
  s1 <- exp(M[i, "theta_1"]); s2 <- exp(M[i, "theta_2"])
  list(b = arr(M[i, "b_x"]), Intercept = M[i, "b_Intercept"] +
         mean(d$x) * M[i, "b_x"], sigma = M[i, "sigma"],
       sd_2 = arr(s1), z_2 = zmat(M, i, "mmg1g2", s1),
       sd_1 = arr(s2), z_1 = zmat(M, i, "g1", s2))
})
bC <- bfix(yp ~ x + (1 | mm(g1, g2)) + (1 | g1), d, inC)
say("C brms r_mmg1g2[2] / r_g1[2] max diff %.3g %.3g",
    mx(as_draws_matrix(bC, variable = "r_mmg1g2[2,Intercept]"), M[, "r_mmg1g2[2,Intercept]"]),
    mx(as_draws_matrix(bC, variable = "r_g1[2,Intercept]"), M[, "r_g1[2,Intercept]"]))
cmp("C mm + (1|g1), g1 = g2 = 2 (no draw)", data.frame(g1 = "2", g2 = "2"), bC, dC, 1)
cmp("C mm + (1|g1), g1 = 2, g2 unset", data.frame(g1 = "2"), bC, dC)
cmp("C mm + (1|g1), nothing set", NULL, bC, dC)
cmp("C mm + (1|g1), g1 = 99, g2 = 2", data.frame(g1 = "99", g2 = "2"), bC, dC)
for (bd in c("wald", "boot")) {
  r <- tryCatch(suppressMessages(conditional_effects(
    fC, "x", resolution = 3, re_formula = NULL, band = bd, boot = 20, seed = 5)),
    error = function(e) e)
  say("C fit, nothing set, band %s: %s", bd, if (inherits(r, "error"))
    paste("ERROR", substr(conditionMessage(r), 1, 160)) else
      paste("width", f4(r$x$upper__ - r$x$lower__)))
}
say("done")
