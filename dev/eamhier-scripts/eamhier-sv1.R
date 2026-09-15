# Lane eamhier: is `sv` recoverable at all, with no hierarchy in the
# way?
#
# Arm C of this lane fits a hierarchical Wiener with across-trial drift
# variability and gets `sv` back at 0.31 and 0.56 against a truth of
# 0.4 (dev/ndt-findings.md, "One thing that DID move"). Two things could
# produce that: the hierarchy, or `sv` itself at these parameter values.
# This script removes the hierarchy and keeps everything else, so that
# the two are separable rather than confounded.
#
# ONE subject, the same 12,000 trials, the same two conditions and the
# same truths. No random effects, so the non-decision bound is the
# global one and it is the right one.
#
# Run:
#   Rscript --vanilla dev/eamhier-scripts/eamhier-sv1.R \
#     <seed> <n> <sv_true> <fit_sv TRUE|FALSE> <outfile>

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5L)
seed <- as.integer(args[[1L]])
n <- as.integer(args[[2L]])
sv_true <- as.numeric(args[[3L]])
fit_sv <- as.logical(args[[4L]])
out <- args[[5L]]

source("dev/eamhier-scripts/eamhier-common.R")
eamhier_libs()
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})

tr <- eamhier_truth
set.seed(seed)
cond <- rep(0:1, each = n / 2L)
d <- ddm_simulate(n, mu = tr$mu0 + tr$mu_cond * cond, bs = tr$bs,
                  ndt = tr$ndt, bias = 0.5, sv = sv_true)
d$cond <- factor(cond, labels = c("a", "b"))

form <- bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1, bias = 0.5)
fam <- if (fit_sv) wiener(variability = "sv") else wiener()

t_start <- Sys.time()
fit <- try(frm(form, family = fam, data = d, se = TRUE), silent = TRUE)
fit_s <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
if (inherits(fit, "try-error")) {
  eamhier_record(out, arm = "S1", seed = seed, n = n,
                 sv_true = sv_true, fit_sv = fit_sv, fit_s = fit_s,
                 status = "error",
                 msg = gsub("[\t\n]", " ", as.character(fit)))
  quit(save = "no")
}

ci <- suppressWarnings(stats::confint(fit))
b <- unlist(fixef(fit))
dg <- frmtmb::diagnose(fit, quiet = TRUE)
nd <- suppressWarnings(
  stats::predict(fit, newdata = d[1L, , drop = FALSE], dpar = "ndt",
                 type = "response", se.fit = TRUE))

eamhier_record(
  out, arm = "S1", seed = seed, n = n, sv_true = sv_true,
  fit_sv = fit_sv, fit_s = fit_s, status = "ok",
  conv = fit$opt$convergence, max_grad = dg$max_grad,
  pdHess = isTRUE(dg$pdHess), n_bad_se = length(dg$bad_se),
  logLik = as.numeric(stats::logLik(fit)),
  mu0 = unname(b[["mu.(Intercept)"]]),
  mu0_lo = eamhier_ci(ci, "(Intercept)")[1L],
  mu0_hi = eamhier_ci(ci, "(Intercept)")[2L],
  mu_cond = unname(b[["mu.condb"]]),
  mu_cond_lo = eamhier_ci(ci, "condb")[1L],
  mu_cond_hi = eamhier_ci(ci, "condb")[2L],
  lbs = unname(b[["bs.(Intercept)"]]),
  lbs_lo = eamhier_ci(ci, "bs_(Intercept)")[1L],
  lbs_hi = eamhier_ci(ci, "bs_(Intercept)")[2L],
  lsv = if (fit_sv) unname(b[["sv.(Intercept)"]]) else NA_real_,
  lsv_lo = eamhier_ci(ci, "sv_(Intercept)")[1L],
  lsv_hi = eamhier_ci(ci, "sv_(Intercept)")[2L],
  ndt_hat = as.numeric(nd$fit[1L]), ndt_se = as.numeric(nd$se.fit[1L]),
  ndt_true = tr$ndt, min_rt = min(d$rt))
