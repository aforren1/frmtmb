# Validation of emmeans support on nonlinear and multivariate fits
# (lane emm, 2026-09-25). Run from the worktree root:
#
#   FRMTMB_LIB=/opt/rlib/lane-emm Rscript dev/emm-validate.R \
#     > dev/emm-validate-log.txt 2>&1
#
# Sections:
#   1. nlpar and whole-mu means against the same model fitted linearly
#   2. resp = and the stacked rep.meas factor against univariate fits
#   3. the delta-method covariance against a numerical Jacobian
#      (numDeriv) and against marginaleffects
#   4. the delta-method standard error against a parametric bootstrap
# Every number below is printed by this script; the findings file
# quotes this log.
lib <- Sys.getenv("FRMTMB_LIB", "/opt/rlib/lane-emm")
.libPaths(c(lib, "/opt/rlib/base", "/opt/rlib/deps", "/opt/r/lib/R/library"))
suppressMessages({
  library(frmtmb)
  library(emmeans)
})
cat("frmtmb from", find.package("frmtmb"), "\n")
B <- as.integer(Sys.getenv("EMM_BOOT", "400"))

sm <- function(e) as.data.frame(summary(e))
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-12))
line <- function(label, x) cat(sprintf("%-58s %.3e\n", label, x))

set.seed(1)
n <- 200
d <- data.frame(f = factor(sample(c("A", "B", "C"), n, TRUE)),
                x = runif(n))
d$y1 <- 1 + as.numeric(d$f) * 0.5 + 2 * d$x + rnorm(n)
d$y2 <- 0.5 * as.numeric(d$f) + rnorm(n)

## 1. nonlinear against linear --------------------------------------------
cat("\n== 1. nlpar and whole mu against the linear fit (seed 1, n = 200)\n")
fnl <- frm(bf(y1 ~ a + b * x, a ~ f, b ~ 1, nl = TRUE), data = d)
flin <- frm(bf(y1 ~ f + x), data = d)
# the two optimizers stop at slightly different points; this is the
# floor every comparison below sits on
cf_nl <- unname(c(fixef(fnl)[c("a_Intercept", "a_fB", "a_fC",
                               "b_Intercept"), "Estimate"]))
cf_lin <- unname(fixef(flin)[c("Intercept", "fB", "fC", "x"), "Estimate"])
line("coefficients, max relative difference (optimizer floor)",
     rel(cf_nl, cf_lin))
a_nl <- sm(emmeans(fnl, "f", nlpar = "a"))
a_lin <- sm(emmeans(flin, "f", at = list(x = 0)))
line("nlpar = 'a' vs linear at x = 0: estimate", rel(a_nl$emmean,
                                                      a_lin$emmean))
line("nlpar = 'a' vs linear at x = 0: SE", rel(a_nl$SE, a_lin$SE))
m_nl <- sm(emmeans(fnl, "f"))
m_lin <- sm(emmeans(flin, "f"))
line("whole mu (delta route) vs linear: estimate", rel(m_nl$emmean,
                                                       m_lin$emmean))
line("whole mu (delta route) vs linear: SE", rel(m_nl$SE, m_lin$SE))
p_nl <- sm(pairs(emmeans(fnl, "f")))
p_lin <- sm(pairs(emmeans(flin, "f")))
line("pairs() on whole mu vs linear: estimate", rel(p_nl$estimate,
                                                    p_lin$estimate))
line("pairs() on whole mu vs linear: SE", rel(p_nl$SE, p_lin$SE))
e_nl <- sm(emmeans(fnl, "f", epred = TRUE))
line("epred = TRUE vs whole mu, identity link (identity): estimate",
     rel(e_nl$emmean, m_nl$emmean))
line("epred = TRUE vs whole mu, identity link (identity): SE",
     rel(e_nl$SE, m_nl$SE))

## 2. multivariate against univariate ---------------------------------------
cat("\n== 2. multivariate against univariate fits (seed 1, n = 200)\n")
fmv <- frm(mvbf(bf(y1 ~ f), bf(y2 ~ f)), data = d)
fu1 <- frm(bf(y1 ~ f), data = d)
fu2 <- frm(bf(y2 ~ f), data = d)
r1 <- sm(emmeans(fmv, "f", resp = "y1"))
u1 <- sm(emmeans(fu1, "f"))
u2 <- sm(emmeans(fu2, "f"))
line("resp = 'y1' vs univariate y1: estimate", rel(r1$emmean, u1$emmean))
line("resp = 'y1' vs univariate y1: SE", rel(r1$SE, u1$SE))
st <- sm(emmeans(fmv, ~ f | rep.meas))
cat("rep.meas levels:", levels(st$rep.meas), "\n")
line("stacked, rep.meas = y1 vs univariate: estimate",
     rel(st$emmean[st$rep.meas == "y1"], u1$emmean))
line("stacked, rep.meas = y2 vs univariate: estimate",
     rel(st$emmean[st$rep.meas == "y2"], u2$emmean))
line("stacked, rep.meas = y2 vs univariate: SE",
     rel(st$SE[st$rep.meas == "y2"], u2$SE))
# without rescor the two responses' coefficients are independent, so
# the y1 - y2 contrast has the root-sum-of-squares standard error
ct <- sm(pairs(emmeans(fmv, ~ rep.meas | f)))
line("contrast y1 - y2 within f: estimate vs difference",
     rel(ct$estimate, u1$emmean - u2$emmean))
line("contrast y1 - y2 within f: SE vs sqrt(se1^2 + se2^2)",
     rel(ct$SE, sqrt(u1$SE^2 + u2$SE^2)))
avg <- sm(emmeans(fmv, "f"))
line("average over rep.meas vs (y1 + y2) / 2: estimate",
     rel(avg$emmean, (u1$emmean + u2$emmean) / 2))
line("average over rep.meas vs (y1 + y2) / 2: SE",
     rel(avg$SE, sqrt(u1$SE^2 + u2$SE^2) / 2))

## 3. delta method against a numerical Jacobian ------------------------------
cat("\n== 3. delta method vs numDeriv and marginaleffects (seed 11)\n")
set.seed(11)
d3 <- data.frame(f = factor(sample(c("A", "B", "C"), n, TRUE)),
                 x = runif(n, 0, 2))
d3$y3 <- c(A = 2, B = 3, C = 4)[as.character(d3$f)] * exp(-0.8 * d3$x) +
  rnorm(n, sd = 0.3)
d3$cnt <- rpois(n, exp(c(A = 0.5, B = 1, C = 1.3)[as.character(d3$f)] -
                         0.6 * d3$x))
fexp <- frm(bf(y3 ~ a * exp(-b * x), a ~ f, b ~ 1, nl = TRUE), data = d3)
fpois <- frm(bf(cnt ~ a - b * x, a ~ f, b ~ 1, nl = TRUE) + poisson(),
             data = d3)

check_jac <- function(fit, rg, type, label) {
  grid <- rg@grid
  th0 <- frmtmb:::interop_coef_vector(fit)
  pred <- function(th) {
    as.numeric(frm_linpred(frmtmb:::set_coef.frmtmb_fit(fit, th),
                           newdata = grid, type = type))
  }
  J <- numDeriv::jacobian(pred, th0)
  Vnum <- J %*% frmtmb:::interop_vcov(fit) %*% t(J)
  line(paste0(label, ": bhat vs predict"), rel(rg@bhat, pred(th0)))
  line(paste0(label, ": V vs numDeriv J V J'"), rel(rg@V, Vnum))
  mfx <- marginaleffects::predictions(fit, newdata = grid, type = type)
  line(paste0(label, ": SE vs marginaleffects"),
       rel(sqrt(diag(rg@V)), mfx$std.error))
}
rg_exp <- ref_grid(fexp, at = list(x = c(0.2, 0.8, 1.6)))
check_jac(fexp, rg_exp, "link", "exp decay, whole mu (6 x 3 grid)")
rg_pl <- ref_grid(fpois, at = list(x = c(0.2, 0.8, 1.6)))
check_jac(fpois, rg_pl, "link", "poisson nl, whole mu, link scale")
rg_pe <- ref_grid(fpois, at = list(x = c(0.2, 0.8, 1.6)), epred = TRUE)
check_jac(fpois, rg_pe, "response", "poisson nl, epred, response scale")
# type = "response" on the link-scale grid applies exp() once, so the
# back-transformed means equal the epred means at each grid point
back <- summary(regrid(rg_pl))
line("poisson nl: regrid(link) vs epred grid: estimate",
     rel(back[[attr(back, "estName")]], rg_pe@bhat))
line("poisson nl: regrid(link) vs epred grid: SE",
     rel(back$SE, sqrt(diag(rg_pe@V))))

## 4. parametric bootstrap --------------------------------------------------
cat("\n== 4. parametric bootstrap, B =", B, "(seed 20260925)\n")
boot_est <- function(fit, sims, dat, resp, fun) {
  out <- vector("list", ncol(sims))
  for (k in seq_len(ncol(sims))) {
    dk <- dat
    dk[[resp]] <- sims[[k]]
    fk <- tryCatch(suppressWarnings(update(fit, data = dk)),
                   error = function(e) NULL)
    out[[k]] <- if (is.null(fk)) NULL else fun(fk)
  }
  do.call(rbind, out)
}
set.seed(20260925)
s_exp <- simulate(fexp, nsim = B, seed = 20260925)
f_exp <- function(ft) {
  c(sm(emmeans(ft, "f"))$emmean, sm(pairs(emmeans(ft, "f")))$estimate)
}
t0 <- proc.time()[["elapsed"]]
bt <- boot_est(fexp, s_exp, d3, "y3", f_exp)
cat("exp decay: refits kept", nrow(bt), "of", B, "in",
    round(proc.time()[["elapsed"]] - t0), "s\n")
se_delta <- c(sm(emmeans(fexp, "f"))$SE, sm(pairs(emmeans(fexp, "f")))$SE)
ratio <- se_delta / apply(bt, 2, sd)
cat("exp decay, whole mu over f then pairs: SE_delta / SD_boot\n")
print(round(ratio, 3))

s_pois <- simulate(fpois, nsim = B, seed = 20260926)
f_pois <- function(ft) {
  c(sm(emmeans(ft, "f", epred = TRUE))$emmean,
    sm(pairs(emmeans(ft, "f", epred = TRUE)))$estimate)
}
t0 <- proc.time()[["elapsed"]]
bp <- boot_est(fpois, s_pois, d3, "cnt", f_pois)
cat("poisson nl: refits kept", nrow(bp), "of", B, "in",
    round(proc.time()[["elapsed"]] - t0), "s\n")
se_delta <- c(sm(emmeans(fpois, "f", epred = TRUE))$SE,
              sm(pairs(emmeans(fpois, "f", epred = TRUE)))$SE)
ratio <- se_delta / apply(bp, 2, sd)
cat("poisson nl, epred over f then pairs: SE_delta / SD_boot\n")
print(round(ratio, 3))
cat("an SD from", B, "draws has relative standard error about",
    round(1 / sqrt(2 * (B - 1)), 3), "\n")
