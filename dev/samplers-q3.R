# Question 3: SparseNUTS (SNUTS) against frm_sample()'s tmbstan route on
# the same frmtmb objective.
#
#   Rscript dev/samplers-q3.R <id> [outdir] [--laplace]
#
# <id> is a design of samplers_q3() in dev/samplers-designs.R.
# Arm "tmbstan": frm_sample(fit) as a user calls it (brms default priors,
#   non-centered blocks where frm_sample() plans them), 4 chains x
#   1000/1000, seed 1, cores 1.
# Arm "snuts": SparseNUTS::sample_snuts() on the objective frm_sample()
#   hands to tmbstan when reparameterize = FALSE (the same priors, the
#   centered blocks: SNUTS decorrelates with the joint precision, which
#   is sparse in the centered coordinates). The objective is captured by
#   a tracer on tmbstan::tmbstan() that stops the call. 4 chains, 1000
#   draws, SNUTS's own warmup default, seed 1, cores 1.
# --laplace: both arms with laplace = TRUE instead.
# Both arms' draws go through frm_sample()'s own map to brms's names
# (frmtmb.sample:::sample_finish()), so they are compared column by
# column.
args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
outdir <- if (length(args) >= 2 && !startsWith(args[2], "--")) args[2] else
  "dev/samplers-q3-out"
lap <- "--laplace" %in% args
tag <- if (lap) paste0(id, "-laplace") else id
source("dev/samplers-designs.R")
samplers_lib()
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
eq <- samplers_eq_env()
suppressPackageStartupMessages({
  library(frmtmb.sample)
  library(SparseNUTS)
  library(posterior)
})
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "| frmtmb.sample", as.character(packageVersion("frmtmb.sample")),
    "| SparseNUTS", as.character(packageVersion("SparseNUTS")),
    "| StanEstimators", as.character(packageVersion("StanEstimators")),
    "| tmbstan", as.character(packageVersion("tmbstan")), "\n")
el <- function() proc.time()[["elapsed"]]
rec <- list(id = id, laplace = lap)
spec <- samplers_design(id, eq)
rec$nobs <- nrow(spec$data)
cat("design", tag, ":", spec$bf, spec$family, "n =", rec$nobs, "\n")

t0 <- el()
fit <- withCallingHandlers(
  tryCatch(eq$eq_frm_fit(spec), error = function(e) e),
  warning = function(w) {
    rec$frm_warn <<- c(rec$frm_warn, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
rec$t_frm <- el() - t0
if (inherits(fit, "error")) {
  rec$frm_error <- conditionMessage(fit)
  cat("frm() error:", rec$frm_error, "\n")
  saveRDS(rec, file.path(outdir, paste0(tag, ".rds")))
  quit(save = "no")
}
rec$n_par <- length(fit$obj$env$last.par.best)
rec$n_random <- length(fit$obj$env$random)

draws_arr <- function(m, nchain) {
  n <- nrow(m) / nchain
  a <- array(NA_real_, c(n, nchain, ncol(m)),
             dimnames = list(NULL, NULL, colnames(m)))
  for (ch in seq_len(nchain)) a[, ch, ] <- m[(ch - 1) * n + seq_len(n), ]
  posterior::as_draws_array(a)
}
summ <- function(da) {
  as.data.frame(posterior::summarise_draws(
    da, "mean", "sd", "mcse_mean", "rhat", "ess_bulk", "ess_tail"))
}

# ---- arm tmbstan: frm_sample() as a user calls it
t0 <- el()
res <- withCallingHandlers(
  tryCatch(frm_sample(fit, chains = 4, iter = 2000, warmup = 1000,
                      seed = 1, cores = 1, refresh = 0, laplace = lap),
           error = function(e) e),
  warning = function(w) {
    rec$tmb_warn <<- c(rec$tmb_warn, conditionMessage(w))
    invokeRestart("muffleWarning")
  },
  message = function(m) {
    rec$tmb_msg <<- c(rec$tmb_msg, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
rec$t_tmb_total <- el() - t0
if (inherits(res, "error")) {
  rec$tmb_error <- conditionMessage(res)
} else {
  sf <- res$stanfit
  rec$t_tmb_sampling <- sum(rstan::get_elapsed_time(sf))
  sp <- rstan::get_sampler_params(sf, inc_warmup = FALSE)
  rec$tmb_div <- sum(vapply(sp, function(x) sum(x[, "divergent__"]), 0))
  rec$tmb_maxtd <- max(vapply(sp, function(x) max(x[, "treedepth__"]), 0))
  rec$tmb_leap <- sum(vapply(sp, function(x) sum(x[, "n_leapfrog__"]), 0))
  rec$tmb_ncp <- res$reparam$labels
  rec$sum_tmb <- summ(draws_arr(res$draws, 4))
}
cat("tmbstan total", rec$t_tmb_total, "sampling", rec$t_tmb_sampling,
    "error", rec$tmb_error %||% "-", "\n")

# ---- capture the objective frm_sample() hands to tmbstan, centered
cap <- new.env()
trace("tmbstan", where = asNamespace("tmbstan"), print = FALSE,
      tracer = bquote({
        assign("obj", obj, envir = .(cap))
        assign("dots", list(...), envir = .(cap))
        stop("samplers: objective captured")
      }))
t0 <- el()
got <- tryCatch(suppressMessages(suppressWarnings(
  frm_sample(fit, chains = 4, iter = 2000, seed = 1, cores = 1,
             refresh = 0, laplace = lap, reparameterize = FALSE))),
  error = function(e) conditionMessage(e))
rec$t_capture <- el() - t0
untrace("tmbstan", where = asNamespace("tmbstan"))
obj <- cap$obj
if (is.null(obj)) {
  rec$capture_error <- got
} else {
  rec$obj_npar <- length(obj$env$last.par.best)
  rec$obj_nouter <- length(obj$par)
  rec$obj_same_as_fit <- identical(obj, fit$obj)
  pn <- names(obj$env$last.par.best)
  t0 <- el()
  sn <- withCallingHandlers(
    tryCatch(sample_snuts(obj, num_samples = 1000, chains = 4, cores = 1,
                          seed = 1, refresh = 0, laplace = lap,
                          model_name = id, print = FALSE),
             error = function(e) e),
    warning = function(w) {
      rec$sn_warn <<- c(rec$sn_warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      rec$sn_msg <<- c(rec$sn_msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  rec$t_sn_total <- el() - t0
  if (inherits(sn, "error")) {
    rec$sn_error <- conditionMessage(sn)
  } else {
    rec$sn_metric <- sn$metric
    rec$sn_warmup <- sn$warmup
    rec$sn_timing <- unlist(sn$timing)
    rec$t_sn_sampling <- sum(sn$timing$time.total)
    spl <- sn$sampler_params
    post <- function(x) x[(nrow(x) - 999):nrow(x), , drop = FALSE]
    rec$sn_div <- sum(vapply(spl, function(x) sum(post(x)[, "divergent__"]),
                             0))
    rec$sn_maxtd <- max(vapply(spl, function(x)
      max(post(x)[, "treedepth__"]), 0))
    rec$sn_leap <- sum(vapply(spl, function(x)
      sum(post(x)[, "n_leapfrog__"]), 0))
    s <- sn$samples
    keep <- -(seq_len(sn$warmup))
    a <- s[keep, , -dim(s)[3], drop = FALSE]
    rec$sn_raw_names <- dimnames(a)[[3]]
    # SNUTS names follow obj$env$last.par.best (laplace: obj$par)
    raw <- if (lap) names(obj$par) else pn
    dimnames(a)[[3]] <- raw
    if (lap) dimnames(a)[[3]] <- raw
    m <- tryCatch(frmtmb.sample:::sample_finish(
      a, fit, list(laplace = lap)), error = function(e) e)
    if (inherits(m, "error")) {
      rec$map_error <- conditionMessage(m)
    } else {
      rec$sum_sn <- summ(draws_arr(m, 4))
    }
  }
}
cat("snuts total", rec$t_sn_total, "sampling", rec$t_sn_sampling,
    "metric", rec$sn_metric %||% "-", "error",
    rec$sn_error %||% rec$capture_error %||% rec$map_error %||% "-", "\n")

# ---- efficiency and agreement
is_re <- function(v) grepl("^(r_|z_|zgp_|s_)", v)
eff <- function(s, t) {
  if (is.null(s)) return(c(NA, NA, NA))
  e <- s$ess_bulk[is.finite(s$ess_bulk)]
  eh <- s$ess_bulk[!is_re(s$variable) & is.finite(s$ess_bulk)]
  c(min(e), min(eh), max(s$rhat, na.rm = TRUE))
}
rec$eff_tmb <- eff(rec$sum_tmb)
rec$eff_sn <- eff(rec$sum_sn)
if (!is.null(rec$sum_tmb) && !is.null(rec$sum_sn)) {
  a <- rec$sum_tmb
  b <- rec$sum_sn[match(a$variable, rec$sum_sn$variable), ]
  ok <- !is.na(b$variable)
  z <- (a$mean[ok] - b$mean[ok]) / sqrt(a$mcse_mean[ok]^2 +
                                          b$mcse_mean[ok]^2)
  names(z) <- a$variable[ok]
  z <- z[is.finite(z)]
  rec$z <- z
  rec$sd_ratio <- stats::setNames(a$sd[ok] / b$sd[ok], a$variable[ok])
  rec$nvar <- length(z)
  rec$only_tmb <- setdiff(a$variable, rec$sum_sn$variable)
  rec$only_sn <- setdiff(rec$sum_sn$variable, a$variable)
}
saveRDS(rec, file.path(outdir, paste0(tag, ".rds")))
num <- function(x) if (is.null(x) || !length(x)) NA_real_ else
  round(as.numeric(x)[1], 3)
hz <- if (!is.null(rec$z)) rec$z[!is_re(names(rec$z))]
row <- data.frame(
  id = tag, n = rec$nobs, n_par = rec$n_par, n_random = rec$n_random,
  t_frm = num(rec$t_frm),
  tmb_total = num(rec$t_tmb_total), tmb_sampling = num(rec$t_tmb_sampling),
  tmb_miness = num(rec$eff_tmb[1]), tmb_miness_hyp = num(rec$eff_tmb[2]),
  tmb_maxrhat = num(rec$eff_tmb[3]), tmb_div = num(rec$tmb_div),
  tmb_maxtd = num(rec$tmb_maxtd), tmb_leap = num(rec$tmb_leap),
  tmb_ess_s = num(rec$eff_tmb[1] / rec$t_tmb_total),
  sn_metric = rec$sn_metric %||% NA, sn_warmup = num(rec$sn_warmup),
  sn_total = num(rec$t_sn_total), sn_sampling = num(rec$t_sn_sampling),
  sn_miness = num(rec$eff_sn[1]), sn_miness_hyp = num(rec$eff_sn[2]),
  sn_maxrhat = num(rec$eff_sn[3]), sn_div = num(rec$sn_div),
  sn_maxtd = num(rec$sn_maxtd), sn_leap = num(rec$sn_leap),
  sn_ess_s = num(rec$eff_sn[1] / rec$t_sn_total),
  nvar = num(rec$nvar), maxz = num(if (length(rec$z)) max(abs(rec$z))),
  maxz_hyp = num(if (length(hz)) max(abs(hz))),
  crit = num(if (length(rec$z)) stats::qnorm(1 - 0.025 / length(rec$z))),
  worst = if (length(rec$z)) names(rec$z)[which.max(abs(rec$z))] else NA,
  sdratio_range = if (length(rec$sd_ratio)) paste(
    round(range(rec$sd_ratio, na.rm = TRUE), 3), collapse = "..") else NA,
  error = substr(paste(c(rec$tmb_error, rec$capture_error, rec$sn_error,
                         rec$map_error), collapse = " | "), 1, 300))
utils::write.table(row, file.path(outdir, paste0(tag, ".tsv")), sep = "\t",
                   row.names = FALSE, quote = TRUE)
print(t(row))
