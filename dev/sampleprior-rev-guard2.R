# Reviewer of lane sampleprior, fix round 2: support_guard() and
# natural_log_jacobian() on doubles and on the tape.
.libPaths(c("C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
core <- asNamespace("frmtmb")
cat("lib:", find.package("frmtmb"), "\n")
pbl <- core$prior_base_logdens
nup <- core$nan_unless_pos

## 1. support_guard(): taken out of prior_base_logdens() by a density
## whose only v-dependence is the guard: uniform(0, 1) at x is
## -log(1) + g(x) + g(1 - x); the reviewer tapes it and reads g
u01 <- core$parse_prior_dist("uniform(0, 1)")
sg_tape <- RTMB::MakeTape(function(x) pbl(x, u01), 0.5)
sg <- function(v) pbl(v, u01)
cat("\n-- support_guard via uniform(0, 1), lower edge --\n")
for (v in c(-1e-3, -1e-200, -1e-290, -1e-299, -1e-301, -5e-324, 0, 1e-310,
            1e-300, 1e-200)) {
  cat(sprintf("x = %-9.3g doubles %-4s tape %-4s tape grad %s\n", v,
              format(sg(v)), format(sg_tape(v)),
              format(sg_tape$jacobian(v))))
}
cat("upper edge: x = 1:", sg(1), sg_tape(1), " x = 1 + 2.2e-16:",
    sg(1 + 2.2e-16), sg_tape(1 + 2.2e-16), "\n")
# Hessian at the edge (rtmb-pitfalls item 22: NaN for x <= 1e-155)
J <- RTMB::MakeTape(function(x) sg_tape$jacobian(x), 0.5)
cat("tape Hessian at x = 0:", J(0), " at 0.5:", J(0.5), " at 1e-150:", J(1e-150), "\n")
# the guard alone, against nan_unless_pos(), for |v| > 1e-284
gd <- function(v) nup(v + 1e-300)
vs <- c(-10^seq(-284, 300, by = 0.5), 10^seq(-284, 300, by = 0.5))
vs <- vs[abs(vs) > 1e-284]
eqd <- identical(gd(vs), nup(vs))
gt <- RTMB::MakeTape(function(v) nup(v + 1e-300), 1)
nt <- RTMB::MakeTape(function(v) nup(v), 1)
vt <- vs[abs(vs) < 9e306]
eqt <- identical(vapply(vt, function(v) gt(v), 1), vapply(vt, function(v) nt(v), 1))
eqg <- identical(vapply(vt, function(v) gt$jacobian(v), 1),
                 vapply(vt, function(v) nt$jacobian(v), 1))
cat(sprintf("guard == nan_unless_pos for %d values |v| > 1e-284: doubles %s, tape %s, tape gradient %s\n",
            length(vs), eqd, eqt, eqg))
# is the + 1e-300 folded away on the tape? if it were, the edge would be NaN
cat("tape of the guard alone at v = 0:", gt(0), " (folded would be NaN)\n")

## every guarded density at its edge(s) and just past, on the tape
cat("\n-- each guarded density on the tape --\n")
dd <- list("gamma(1, 1)" = c(0, -1e-12), "gamma(2, 1)" = c(0, -1e-12),
           "exponential(2)" = c(0, -1e-12),
           "inv_gamma(0.4, 0.3)" = c(0, -1e-12),
           "lognormal(0, 1)" = c(0, -1e-12),
           "beta(1, 1)" = c(0, -1e-12, 1, 1 + 1e-12),
           "beta(2, 5)" = c(0, -1e-12, 1, 1 + 1e-12),
           "uniform(0.2, 0.7)" = c(0.2, 0.2 - 1e-12, 0.7, 0.7 + 1e-12))
for (nm in names(dd)) {
  dist <- core$parse_prior_dist(nm)
  tp <- RTMB::MakeTape(function(x) pbl(x, dist), 0.5)
  v <- dd[[nm]]
  cat(sprintf("%-20s %s\n", nm, paste(sprintf("x=%-.12g: %s (grad %s)", v,
    format(sapply(v, function(z) tp(z)), digits = 6),
    format(sapply(v, function(z) tp$jacobian(z)), digits = 4)), collapse = " | ")))
}

## 2. natural_log_jacobian(): logit against plogis(log.p = TRUE) and a
## 256-bit reference, value and gradient, doubles and tape
cat("\n-- natural_log_jacobian, logit --\n")
lk <- frmtmb:::get_link("logit")
nlj <- core$natural_log_jacobian
eta <- seq(-50, 50, by = 0.25)
ref <- vapply(eta, function(e) {
  E <- Rmpfr::mpfr(e, 256)
  as.numeric(-log1p(exp(-E)) - log1p(exp(E)))
}, 1)
refg <- vapply(eta, function(e) {
  E <- Rmpfr::mpfr(e, 256)
  as.numeric(1 - 2 / (1 + exp(-E)))
}, 1)
pl <- stats::plogis(eta, log.p = TRUE) + stats::plogis(-eta, log.p = TRUE)
dbl <- vapply(eta, function(e) nlj(e, lk), 1)
tp <- RTMB::MakeTape(function(e) nlj(e, lk), 0)
tv <- vapply(eta, function(e) tp(e), 1)
tg <- vapply(eta, function(e) tp$jacobian(e), 1)
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))
cat(sprintf("value: doubles vs 256-bit max rel %.2g; tape vs 256-bit %.2g; plogis(log.p) vs 256-bit %.2g\n",
            rel(dbl, ref), rel(tv, ref), rel(pl, ref)))
cat(sprintf("gradient: tape vs 256-bit (1 - 2p) max abs %.2g (max |ref| 1)\n",
            max(abs(tg - refg))))
old <- log(abs(lk$mu_eta(eta)))
cat(sprintf("0.71.0's log(mu_eta) is -Inf at %d of %d points (|eta| >= %.2f)\n",
            sum(!is.finite(old)), length(eta), min(abs(eta[!is.finite(old)]))))

## frmtmb.eam's scaled_logit (wiener ndt): top * plogis(eta)
set.seed(20261115)
dw <- frmtmb.eam::ddm_simulate(50, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
fit <- frm(rt | dec(upper) ~ 1, data = dw, family = frmtmb.eam::wiener(),
           dry_run = "objective")
lp <- Filter(function(l) identical(l$dpar, "ndt"), fit$frame$linpreds)[[1]]
sl <- lp$link
top <- sl$linkinv(Inf)
cat("\n-- scaled_logit:", sl$name, " top", top, " min rt", min(dw$rt),
    " linkinv == top * plogis:", isTRUE(all.equal(sl$linkinv(eta),
                                                 top * plogis(eta), tolerance = 0)), "\n")
refs <- vapply(eta, function(e) {
  E <- Rmpfr::mpfr(e, 256)
  as.numeric(log(Rmpfr::mpfr(top, 256)) - log1p(exp(-E)) - log1p(exp(E)))
}, 1)
tps <- RTMB::MakeTape(function(e) nlj(e, sl), 0)
cat(sprintf("value: doubles vs 256-bit max rel %.2g; tape %.2g; gradient tape vs 1 - 2p max abs %.2g\n",
            rel(vapply(eta, function(e) nlj(e, sl), 1), refs),
            rel(vapply(eta, function(e) tps(e), 1), refs),
            max(abs(vapply(eta, function(e) tps$jacobian(e), 1) - refg))))

## every other link: natural_log_jacobian() is log|mu_eta|, bitwise
cat("\n-- other links --\n")
for (nm in c("log", "identity", "softplus", "squareplus", "logm1", "log1p",
             "inverse", "sqrt")) {
  l <- tryCatch(frmtmb:::get_link(nm), error = function(e) NULL)
  if (is.null(l)) { cat(nm, ": no such link\n"); next }
  e <- seq(-5, 5, by = 0.5)
  a <- suppressWarnings(vapply(e, function(z) nlj(z, l), 1))
  b <- suppressWarnings(log(abs(l$mu_eta(e))))
  cat(sprintf("%-11s identical to log|mu_eta|: %s\n", nm, identical(a, b)))
}
