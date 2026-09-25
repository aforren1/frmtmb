# Lane arcov: brms's cov = FALSE ARMA, validated three ways.
#
#  1. The log-likelihood at a shared parameter point against an R
#     transliteration of the model block brms 2.23.0 generates
#     (brms::make_stancode() read by hand), fed with brms's OWN Stan data
#     (brms::make_standata(): Y in brms's (gr, time) order and J_lag).
#     Nothing on the reference side calls frmtmb's recursion or its
#     row ordering. Points: the ML optimum and a random perturbation of
#     it (seeded), so agreement is not an artifact of a zero gradient.
#  2. The ML estimates of a single long series against the conditional
#     sum of squares, stats::arima(method = "CSS"), for AR(1) and MA(1).
#  3. Fit time at N = 2000 in 100 groups.
#
# Run: Rscript dev/arcov-validate.R  (from the worktree root; prints the
# table recorded in dev/arcov-findings.md)
.libPaths(c("/opt/rlib/lane-arcov", "/opt/rlib/base", "/opt/rlib/deps",
            "/opt/r/lib/R/library"))
# brms first, so that the bare bf() below is frmtmb's; brms's own
# formulas for make_standata() are written brms::bf()
suppressMessages({
  library(brms)
  library(frmtmb)
})

# ---- brms's model block, transliterated -------------------------------
# for (n in 1:N) {
#   mu[n] += Err[n, 1:Kma] * ma;
#   err[n] = Y[n] - mu[n];
#   for (i in 1:J_lag[n]) Err[n + 1, i] = err[n + 1 - i];
#   mu[n] += Err[n, 1:Kar] * ar;
# }
brms_arma_mu <- function(mu, Y, J_lag, ar, ma) {
  N <- length(Y)
  Kar <- length(ar)
  Kma <- length(ma)
  max_lag <- max(Kar, Kma)
  Err <- matrix(0, N + 1, max_lag)
  err <- numeric(N)
  for (n in seq_len(N)) {
    if (Kma) mu[n] <- mu[n] + sum(Err[n, seq_len(Kma)] * ma)
    err[n] <- Y[n] - mu[n]
    for (i in seq_len(J_lag[n])) Err[n + 1, i] <- err[n + 1 - i]
    if (Kar) mu[n] <- mu[n] + sum(Err[n, seq_len(Kar)] * ar)
  }
  mu
}

# natural dpar values of a frmtmb fit at an internal parameter list, in
# frmtmb's row order; the parameterization only, no likelihood
dpars_at <- function(f, pl) {
  out <- list()
  for (lp in f$frame$linpreds) {
    co <- pl[[lp$par]][lp$idx]
    eta <- if (ncol(lp$X)) drop(as.matrix(lp$X) %*% co) else
      numeric(f$frame$n_obs)
    if (!is.null(lp$Z)) eta <- eta + drop(as.matrix(lp$Z %*% pl$b))
    out[[lp$resp]][[lp$dpar]] <- lp$link$linkinv(eta)
  }
  out
}

# the internal parameter list at a full parameter vector
par_list <- function(f, full) {
  nm <- names(full)
  lapply(split(unname(full), factor(nm, unique(nm))), identity)
}

# the joint negative log-likelihood frmtmb tapes, at a full vector
# (random effects included, so no Laplace step is involved)
frm_joint_nll <- function(f, full) f$obj$env$f(full, order = 0)

# one case: fit, then compare at the optimum and at a perturbed point
check_case <- function(label, formula, data, family, ref_ll, seed = 1,
                       fam_brms = brms::brmsfamily(family[["family"]]),
                       bform = formula, ...) {
  f <- frm(formula, data = data, family = family, ...)
  sd <- make_standata(bform, data = data, family = fam_brms,
                      internal = TRUE)
  ord <- order(attr(sd, "old_order"))
  full0 <- f$obj$env$last.par.best
  set.seed(seed)
  fullp <- full0 + stats::rnorm(length(full0), 0, 0.05)
  res <- NULL
  for (pt in c("optimum", "perturbed")) {
    full <- if (pt == "optimum") full0 else fullp
    pl <- par_list(f, full)
    ll_ref <- ref_ll(f, sd, ord, pl)
    ll_frm <- -frm_joint_nll(f, full)
    res <- rbind(res, data.frame(
      case = label, point = pt, loglik_frmtmb = ll_frm,
      loglik_brms_recursion = ll_ref,
      rel_diff = abs(ll_frm - ll_ref) / abs(ll_ref)))
  }
  res
}

# ---- data --------------------------------------------------------------
set.seed(20260925)
d <- expand.grid(week = 1:8, subj = factor(1:12))
# ragged, with interior gaps, and shuffled: brms counts lags in rows
# within (gr, time) order, and the data order must not matter
d <- d[-c(3, 17, 18, 40, 41, 42, 77, 90), ]
d$x <- stats::rnorm(nrow(d))
d$y <- 1 + 0.5 * d$x + stats::rnorm(nrow(d))
d$w <- stats::runif(nrow(d), 0.5, 2)
d$cc <- sample(c(0, 0, 0, 1, -1), nrow(d), TRUE)
d$yp <- d$y + 3
d <- d[sample(nrow(d)), ]
rownames(d) <- NULL

coef_ar <- function(f, pl) {
  ac <- f$frame$autocor[[1]]
  th <- pl$thetaac[ac$theta_idx]
  list(ar = th[seq_len(ac$p)], ma = th[ac$p + seq_len(ac$q)])
}

ref_gauss <- function(f, sd, ord, pl) {
  dp <- dpars_at(f, pl)[[1]]
  cf <- coef_ar(f, pl)
  stopifnot(isTRUE(all.equal(as.numeric(sd$Y),
                             f$frame$y[[1]][ord], tolerance = 0)))
  mu <- brms_arma_mu(rep_len(dp$mu, length(ord))[ord], sd$Y, sd$J_lag,
                     cf$ar, cf$ma)
  sum(stats::dnorm(sd$Y, mu, rep_len(dp$sigma, length(ord))[ord],
                   log = TRUE))
}

ref_student <- function(f, sd, ord, pl) {
  dp <- dpars_at(f, pl)[[1]]
  cf <- coef_ar(f, pl)
  mu <- brms_arma_mu(rep_len(dp$mu, length(ord))[ord], sd$Y, sd$J_lag,
                     cf$ar, cf$ma)
  nu <- rep_len(dp$nu, length(ord))[ord]
  sg <- rep_len(dp$sigma, length(ord))[ord]
  sum(stats::dt((sd$Y - mu) / sg, nu, log = TRUE) - log(sg))
}

ref_weights <- function(f, sd, ord, pl) {
  dp <- dpars_at(f, pl)[[1]]
  cf <- coef_ar(f, pl)
  mu <- brms_arma_mu(rep_len(dp$mu, length(ord))[ord], sd$Y, sd$J_lag,
                     cf$ar, cf$ma)
  sum(sd$weights * stats::dnorm(sd$Y, mu, dp$sigma, log = TRUE))
}

ref_cens <- function(f, sd, ord, pl) {
  dp <- dpars_at(f, pl)[[1]]
  cf <- coef_ar(f, pl)
  mu <- brms_arma_mu(rep_len(dp$mu, length(ord))[ord], sd$Y, sd$J_lag,
                     cf$ar, cf$ma)
  s <- dp$sigma
  cens <- sd$cens
  sum(ifelse(cens == 0, stats::dnorm(sd$Y, mu, s, log = TRUE),
             ifelse(cens == 1,
                    stats::pnorm(sd$Y, mu, s, lower.tail = FALSE,
                                 log.p = TRUE),
                    stats::pnorm(sd$Y, mu, s, log.p = TRUE))))
}

ref_trunc <- function(f, sd, ord, pl) {
  dp <- dpars_at(f, pl)[[1]]
  cf <- coef_ar(f, pl)
  mu <- brms_arma_mu(rep_len(dp$mu, length(ord))[ord], sd$Y, sd$J_lag,
                     cf$ar, cf$ma)
  s <- dp$sigma
  sum(stats::dnorm(sd$Y, mu, s, log = TRUE) -
        stats::pnorm(sd$lb, mu, s, lower.tail = FALSE, log.p = TRUE))
}

# a random intercept: the row density at the conditional mu, plus the
# block density, which is dnorm(b, 0, exp(theta)) for a 1-d us block
ref_re <- function(f, sd, ord, pl) {
  ref_gauss(f, sd, ord, pl) +
    sum(stats::dnorm(pl$b, 0, exp(pl$theta), log = TRUE))
}

ref_rescor <- function(f, sd, ord, pl) {
  dpa <- dpars_at(f, pl)
  acs <- f$frame$autocor
  rs <- names(f$spec$responses)
  Z <- NULL
  lsig <- 0
  for (r in rs) {
    th <- pl$thetaac[acs[[r]]$theta_idx]
    Yr <- sd[[paste0("Y_", r)]]
    mu <- brms_arma_mu(dpa[[r]]$mu[ord], Yr, sd[[paste0("J_lag_", r)]],
                       th[seq_len(acs[[r]]$p)], numeric(0))
    sg <- rep_len(dpa[[r]]$sigma, length(ord))[ord]
    Z <- cbind(Z, (Yr - mu) / sg)
    lsig <- lsig + sum(log(sg))
  }
  C <- frmtmb:::us_chol_cor(pl$thetar, 2)
  sum(mvtnorm::dmvnorm(Z, sigma = C, log = TRUE)) - lsig
}

out <- rbind(
  check_case("gaussian ar(1)", y ~ x + ar(week, subj), d, gaussian(),
             ref_gauss),
  check_case("gaussian ma(1)", y ~ x + ma(week, subj), d, gaussian(),
             ref_gauss),
  check_case("gaussian arma(1,1)", y ~ x + arma(week, subj), d,
             gaussian(), ref_gauss),
  check_case("gaussian arma(2,2)", y ~ x + arma(week, subj, p = 2, q = 2),
             d, gaussian(), ref_gauss),
  check_case("gaussian ar(3)", y ~ x + ar(week, subj, p = 3), d,
             gaussian(), ref_gauss),
  check_case("gaussian ma(2), no time", y ~ x + ma(gr = subj, q = 2), d,
             gaussian(), ref_gauss),
  check_case("gaussian arma(1,1), sigma ~ x",
             bf(y ~ x + arma(week, subj), sigma ~ x), d, gaussian(),
             ref_gauss,
             bform = brms::bf(y ~ x + arma(week, subj), sigma ~ x)),
  check_case("gaussian(log) ar(1)", yp ~ x + ar(week, subj), d,
             gaussian(link = "log"), ref_gauss,
             fam_brms = brms::brmsfamily("gaussian", link = "log")),
  check_case("student arma(1,1)", y ~ x + arma(week, subj), d, student(),
             ref_student),
  check_case("student arma(1,1), nu ~ x",
             bf(y ~ x + arma(week, subj), nu ~ x), d, student(),
             ref_student,
             bform = brms::bf(y ~ x + arma(week, subj), nu ~ x)),
  check_case("weights ma(1)", y | weights(w) ~ x + ma(week, subj), d,
             gaussian(), ref_weights),
  check_case("cens arma(1,1)", y | cens(cc) ~ x + arma(week, subj), d,
             gaussian(), ref_cens),
  check_case("trunc ma(1)", yp | trunc(lb = 0) ~ x + ma(week, subj), d,
             gaussian(), ref_trunc),
  check_case("(1 | subj) + arma(1,1)",
             y ~ x + (1 | subj) + arma(week, subj), d, gaussian(), ref_re)
)
# mi() on the response: brms's Yl, the observed-or-imputed value, is
# what the residual is taken against, and the imputed values are
# parameters of the joint density evaluated here
ref_mi <- function(f, sd, ord, pl) {
  dp <- dpars_at(f, pl)[[1]]
  cf <- coef_ar(f, pl)
  mm <- f$frame$mi_map[[1]]
  yl <- f$frame$y[[1]]
  yl[mm$rows] <- pl$miss[mm$idx]
  Yl <- yl[ord]
  mu <- brms_arma_mu(rep_len(dp$mu, length(ord))[ord], Yl, sd$J_lag,
                     cf$ar, cf$ma)
  sum(stats::dnorm(Yl, mu, dp$sigma, log = TRUE))
}
dm <- d
dm$y[c(5, 20, 33)] <- NA
out <- rbind(out, check_case(
  "mi() response, arma(1,1)", bf(y | mi() ~ x + arma(week, subj)), dm,
  gaussian(), ref_mi, bform = brms::bf(y | mi() ~ x + arma(week, subj))))
d2 <- d
d2$y2 <- d2$y + stats::rnorm(nrow(d2))
out <- rbind(out, check_case(
  "rescor, ar(1) on both",
  bf(y ~ x + ar(week, subj)) + bf(y2 ~ x + ar(week, subj)) +
    set_rescor(TRUE), d2, gaussian(), ref_rescor,
  bform = brms::bf(y ~ x + ar(week, subj)) +
    brms::bf(y2 ~ x + ar(week, subj)) + brms::set_rescor(TRUE)))
old <- options(width = 120, digits = 12)
print(out, row.names = FALSE)
cat("max rel_diff:", format(max(out$rel_diff), digits = 3), "\n")
options(old)

# ---- 2. stats::arima(method = "CSS") on one long series ---------------
# CSS drops the first p rows of an AR(p); brms keeps row 1 with no
# lagged term. Weight 0 on that row removes its density and keeps its
# residual as the lag of row 2, which is exactly CSS. A pure MA(q)
# conditions on nothing in either (e_0 = 0), so no weight is needed.
set.seed(7)
n <- 2000
s <- data.frame(t = seq_len(n),
                y_ar = 2 + as.numeric(stats::arima.sim(list(ar = 0.6), n)),
                y_ma = 1 + as.numeric(stats::arima.sim(list(ma = 0.5), n)))
s$w <- c(0, rep(1, n - 1))
tight <- list(reltol = 1e-15, maxit = 5000)
# both optimizers are run tight: nlminb's default relative tolerance
# stops this N = 2000 fit at a gradient of about 2e-3
fctl <- frmtmb_control(optCtrl = list(iter.max = 1000, eval.max = 1000,
                                      rel.tol = 1e-14))
fa <- frm(y_ar | weights(w) ~ 1 + ar(t), data = s, family = gaussian(),
          control = fctl)
ca <- stats::arima(s$y_ar, order = c(1, 0, 0), method = "CSS",
                   optim.control = tight)
fm <- frm(y_ma ~ 1 + ma(t), data = s, family = gaussian(), control = fctl)
cm <- stats::arima(s$y_ma, order = c(0, 0, 1), method = "CSS",
                   optim.control = tight)
css <- rbind(
  data.frame(model = "AR(1)", par = c("mean", "ar1", "sigma2"),
             frmtmb = c(fa$estimates$beta, fa$estimates$thetaac,
                        exp(2 * fa$estimates$betad)),
             arima_css = c(ca$coef[["intercept"]], ca$coef[["ar1"]],
                           ca$sigma2)),
  data.frame(model = "MA(1)", par = c("mean", "ma1", "sigma2"),
             frmtmb = c(fm$estimates$beta, fm$estimates$thetaac,
                        exp(2 * fm$estimates$betad)),
             arima_css = c(cm$coef[["intercept"]], cm$coef[["ma1"]],
                           cm$sigma2)))
css$rel_diff <- abs(css$frmtmb - css$arima_css) / abs(css$arima_css)
old <- options(width = 120, digits = 12)
print(css, row.names = FALSE)
# The objectives themselves, at arima's estimates: frmtmb's
# log-likelihood there against the gaussian log-likelihood of arima's
# own CSS residuals (its first ncond residuals are the conditioned-on
# rows, which carry weight 0 above). This is where the two optimizers
# cannot disagree.
at_css <- function(f, a, par_nm, drop) {
  pr <- f$obj$par
  pr[names(pr) == "beta"] <- a$coef[["intercept"]]
  pr[names(pr) == "betad"] <- 0.5 * log(a$sigma2)
  pr[names(pr) == "thetaac"] <- a$coef[[par_nm]]
  e <- stats::residuals(a)
  if (drop) e <- e[-1L]
  c(frmtmb = -f$obj$fn(pr),
    arima_resid = sum(stats::dnorm(e, 0, sqrt(a$sigma2), log = TRUE)))
}
obj <- rbind(`AR(1)` = at_css(fa, ca, "ar1", TRUE),
             `MA(1)` = at_css(fm, cm, "ma1", FALSE))
print(cbind(obj, rel_diff = abs(obj[, 1] - obj[, 2]) / abs(obj[, 2])))
options(old)

# ---- 3. fit time, N = 2000 in 100 groups of 20 -------------------------
set.seed(11)
tg <- expand.grid(week = 1:20, subj = factor(1:100))
tg$x <- stats::rnorm(nrow(tg))
tg$y <- 1 + 0.5 * tg$x + as.numeric(replicate(100, stats::arima.sim(
  list(ar = 0.5, ma = 0.3), 20)))
time_fit <- function(fo) {
  best <- Inf
  for (k in 1:3) {
    t0 <- proc.time()[["elapsed"]]
    frm(fo, data = tg, family = gaussian())
    best <- min(best, proc.time()[["elapsed"]] - t0)
  }
  best
}
tm <- data.frame(
  model = c("control: y ~ x (no autocor)", "ar(week, subj)",
            "ma(week, subj)", "arma(week, subj)",
            "arma(week, subj, cov = TRUE)"),
  min_of_3_s = c(time_fit(y ~ x), time_fit(y ~ x + ar(week, subj)),
                 time_fit(y ~ x + ma(week, subj)),
                 time_fit(y ~ x + arma(week, subj)),
                 time_fit(y ~ x + arma(week, subj, cov = TRUE))))
print(tm, row.names = FALSE)

# Per-evaluation cost of the taped objective and gradient, which is what
# an optimizer pays: blocks of at least 1.2 s, the arms interleaved in
# one process, the minimum of 5 rounds. The control is the plain model;
# "control again" is the same object timed a second time and should
# read 1.00. The last arm is the worst shape for the MA recursion, one
# series of 2000 rows, so 2000 positions of one row each.
one <- data.frame(t = seq_len(2000), x = stats::rnorm(2000))
one$y <- 1 + 0.5 * one$x + as.numeric(stats::arima.sim(list(ma = 0.4), 2000))
objs <- list(
  control = frm(y ~ x, data = tg, family = gaussian())$obj,
  ar = frm(y ~ x + ar(week, subj), data = tg, family = gaussian())$obj,
  ma = frm(y ~ x + ma(week, subj), data = tg, family = gaussian())$obj,
  arma = frm(y ~ x + arma(week, subj), data = tg, family = gaussian())$obj,
  arma_cov_true = frm(y ~ x + arma(week, subj, cov = TRUE), data = tg,
                      family = gaussian())$obj,
  ma_one_series = frm(y ~ x + ma(t), data = one, family = gaussian())$obj)
objs <- c(objs[1], list(control_again = objs[[1]]), objs[-1])
per_eval <- function(o, reps) {
  p0 <- o$par
  t0 <- proc.time()[["elapsed"]]
  for (k in seq_len(reps)) {
    o$fn(p0)
    o$gr(p0)
  }
  (proc.time()[["elapsed"]] - t0) / reps
}
reps <- vapply(objs, function(o) {
  r <- 50L
  while (per_eval(o, r) * r < 1.2) r <- r * 2L
  r
}, 1L)
best <- rep(Inf, length(objs))
for (round in 1:5) {
  for (k in seq_along(objs)) {
    best[k] <- min(best[k], per_eval(objs[[k]], reps[[k]]))
  }
}
pe <- data.frame(model = names(objs), us_per_fn_plus_gr = 1e6 * best,
                 ratio_to_control = best / best[1])
print(pe, row.names = FALSE, digits = 4)
