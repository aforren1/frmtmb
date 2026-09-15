# Lane eamhier, item 2.1: ONE replicate, in its own process.
#
# One process per replicate rather than a loop, because a 12,000-row
# hierarchical DDM peaks near 5 GB of process working set and a worker
# that keeps several tapes alive is a machine that swaps. The fit costs
# about 150 s and a fresh R costs about 3 s, so the process is 2 percent
# of the replicate.
#
# Run:
#   Rscript --vanilla dev/eamhier-scripts/eamhier-rep.R \
#     <arm A|B|C> <bound pg|gl> <seed> <ns> <nt> <outfile>
#
# Writes ONE tab-separated key=value line to <outfile>. One file per
# replicate: parallel workers never share a file.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 6L)
arm <- args[[1L]]
bound <- args[[2L]]
seed <- as.integer(args[[3L]])
ns <- as.integer(args[[4L]])
nt <- as.integer(args[[5L]])
out <- args[[6L]]

source("dev/eamhier-scripts/eamhier-common.R")
eamhier_libs()
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})

d <- eamhier_data(seed, arm, ns = ns, nt = nt)
form <- eamhier_form(arm, bound)
peak_ws <- eamhier_peak_ws
fam <- eamhier_family(arm)
t0_true <- attr(d, "ndt_subject")
key <- d[match(levels(d$s), as.character(d$s)), , drop = FALSE]
floors <- as.numeric(tapply(d$rt, d$s, min))[
  match(as.character(key$s), levels(d$s))]

gc(reset = TRUE)
t_start <- Sys.time()
fit <- try(frm(form, family = fam, data = d, se = TRUE), silent = TRUE)
fit_s <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

if (inherits(fit, "try-error")) {
  eamhier_record(out, arm = arm, bound = bound, seed = seed, ns = ns,
                 nt = nt, rows = nrow(d), fit_s = fit_s,
                 status = "error",
                 msg = gsub("[\t\n]", " ", as.character(fit)))
  quit(save = "no")
}

ci <- suppressWarnings(stats::confint(fit))
b <- unlist(fixef(fit))
dg <- frmtmb::diagnose(fit, quiet = TRUE)
vc <- VarCorr(fit)
vsd <- vapply(vc, function(m) sqrt(m[1L, 1L]), numeric(1))

# `theta_k` is the log standard deviation of the k-th variance
# component, in the order VarCorr() lists them. The record carries BOTH
# spellings so that the assumed order is checkable in the file rather
# than in this comment.
th <- function(k) eamhier_ci(ci, paste0("theta_", k))
sd_ci <- function(k) exp(th(k))
th_est <- function(k) {
  j <- which(rownames(ci) == paste0("theta_", k))
  if (!length(j)) NA_real_ else as.numeric(ci[j, 3L])
}

# The population non-decision time in seconds.
#
# The two parameterizations report `ndt` on two different scales and
# only one of them needs the bound put back. Under `ndt_group()` the
# link is a plain logit on a FRACTION of the row's own bound, so the
# population value is that fraction on the mean of the bounds the fit
# used. Under the scalar bound the bound stays INSIDE the link, a
# scaled logit onto (0, ub), so `predict()` already reports a time and
# multiplying by `ub` again would report `time * ub`. Measured, not
# assumed: dev/eamhier-scripts/eamhier-scalar-link.R prints
# `linkinv(eta)`, `predict()` and `ndt_time()` side by side for both,
# and tests/testthat/test-scale.R got this wrong for its ungrouped row.
nd <- suppressWarnings(
  stats::predict(fit, newdata = key[1L, , drop = FALSE], dpar = "ndt",
                 type = "response", re.form = NA, se.fit = TRUE))
bnd <- frmtmb::single_response(fit)[["family"]][["ndt_bound"]]
fl <- if (is.null(bnd[["floors"]])) 1 else mean(bnd[["floors"]])
ndt_frac <- as.numeric(nd$fit[1L])
ndt_frac_se <- as.numeric(nd$se.fit[1L])
ndt_hat <- ndt_frac * fl
ndt_se <- ndt_frac_se * fl

hat <- as.numeric(ndt_time(fit, newdata = key))
err_ms <- 1000 * (hat - t0_true)
margin_ms <- 1000 * (floors - hat)

eamhier_record(
  out,
  arm = arm, bound = bound, seed = seed, ns = ns, nt = nt,
  rows = nrow(d), n_par = length(fit$opt$par), fit_s = fit_s,
  status = "ok",
  conv = fit$opt$convergence, max_grad = dg$max_grad,
  pdHess = isTRUE(dg$pdHess), n_bad_se = length(dg$bad_se),
  n_flat = length(dg$flat),
  n_extreme_theta = length(dg$extreme_theta),
  unbounded = if (is.null(dg$unbounded_dpar)) "none" else
    paste(dg$unbounded_dpar$parameter, collapse = "/"),
  n_sing = length(dg$singular),
  logLik = as.numeric(stats::logLik(fit)),
  # fixed effects, on the scale confint() reports them
  mu0 = unname(b[["mu.(Intercept)"]]),
  mu0_lo = eamhier_ci(ci, "(Intercept)")[1L],
  mu0_hi = eamhier_ci(ci, "(Intercept)")[2L],
  mu_cond = if (identical(arm, "B")) NA_real_ else
    unname(b[["mu.condb"]]),
  mu_cond_lo = eamhier_ci(ci, "condb")[1L],
  mu_cond_hi = eamhier_ci(ci, "condb")[2L],
  lbs = unname(b[["bs.(Intercept)"]]),
  lbs_lo = eamhier_ci(ci, "bs_(Intercept)")[1L],
  lbs_hi = eamhier_ci(ci, "bs_(Intercept)")[2L],
  lsv = if (identical(arm, "C")) unname(b[["sv.(Intercept)"]]) else
    NA_real_,
  lsv_lo = eamhier_ci(ci, "sv_(Intercept)")[1L],
  lsv_hi = eamhier_ci(ci, "sv_(Intercept)")[2L],
  # variance components, as log standard deviations and as standard
  # deviations, in VarCorr()'s order
  vc_names = paste(names(vc), collapse = ";"),
  vc_sd = paste(formatC(vsd, digits = 10, format = "g"),
                collapse = ";"),
  th1 = th_est(1L), th1_lo = th(1L)[1L], th1_hi = th(1L)[2L],
  th2 = th_est(2L), th2_lo = th(2L)[1L], th2_hi = th(2L)[2L],
  th3 = th_est(3L), th3_lo = th(3L)[1L], th3_hi = th(3L)[2L],
  sd1 = exp(th_est(1L)), sd1_lo = sd_ci(1L)[1L],
  sd1_hi = sd_ci(1L)[2L],
  sd2 = exp(th_est(2L)), sd2_lo = sd_ci(2L)[1L],
  sd2_hi = sd_ci(2L)[2L],
  sd3 = exp(th_est(3L)), sd3_lo = sd_ci(3L)[1L],
  sd3_hi = sd_ci(3L)[2L],
  # the non-decision time, which is the quantity item 1.0a's note says
  # to read per subject rather than as a spread
  ndt_frac = ndt_frac, ndt_frac_se = ndt_frac_se,
  ndt_hat = ndt_hat, ndt_se = ndt_se,
  ndt_pop_true = mean(t0_true),
  ndt_rmse_ms = sqrt(mean(err_ms^2)),
  ndt_bias_ms = mean(err_ms),
  ndt_maxabs_ms = max(abs(err_ms)),
  ndt_corr = if (stats::sd(t0_true) > 0) {
    stats::cor(hat, t0_true)
  } else NA_real_,
  ndt_sd_hat_ms = 1000 * stats::sd(hat),
  ndt_sd_true_ms = 1000 * stats::sd(t0_true),
  ndt_margin_min_ms = min(margin_ms),
  n_below_own_floor = sum(margin_ms > 0),
  floor_sd_ms = 1000 * stats::sd(floors),
  floor_min_ms = 1000 * min(floors),
  floor_mean_ms = 1000 * mean(floors),
  # the truths the draw actually used, so that a replicate is judged
  # against its own draw where that is the right comparison
  sd_mu_drawn = stats::sd(attr(d, "mu_subject")),
  sd_bs_drawn = stats::sd(log(attr(d, "bs_subject"))),
  # The R heap is a FLOOR on what a fit costs: RTMB's tape lives in C++
  # memory that gc() cannot see, so the peak process working set rides
  # along beside it. eamhier-memory.R has the stage-by-stage split.
  mem_mb = sum(gc()[, 6L]),
  peak_ws_mb = peak_ws())
