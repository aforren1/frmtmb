# Question 1: brms + stanli against brms + rstan (CmdStan's sampler in
# rstan) and against frm(), on one design.
#
#   Rscript dev/samplers-q1.R <id> [outdir]
#
# <id> is a survey pair (pNNN) or a large design (big_*), see
# dev/samplers-designs.R. The brms side uses brms's DEFAULT priors (plus
# a prior the spec itself sets); frm() is the maximum likelihood fit.
# Both samplers run 4 chains x 1000 warmup / 1000 draws, sequentially,
# seed 1. rstan compiles fresh (no cache), so its time is what a new
# model costs. Writes <outdir>/<id>.rds and <outdir>/<id>.tsv.
args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
outdir <- if (length(args) >= 2) args[2] else "dev/samplers-q1-out"
source("dev/samplers-designs.R")
samplers_lib()
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
eq <- samplers_eq_env()
suppressPackageStartupMessages({
  library(stanli)
  library(posterior)
})
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "| stanli", as.character(packageVersion("stanli")),
    "| rstan", as.character(packageVersion("rstan")),
    "| brms", as.character(packageVersion("brms")), "\n")
el <- function() proc.time()[["elapsed"]]
rec <- list(id = id)
spec <- samplers_design(id, eq)
rec$key <- spec$key
rec$nobs <- nrow(spec$data)
cat("design", id, ":", spec$key, "n =", rec$nobs, "\n")

# ---- frm(): maximum likelihood
t0 <- el()
fit <- tryCatch(eq$eq_frm_fit(spec), error = function(e) e)
rec$t_frm <- el() - t0
if (inherits(fit, "error")) {
  rec$frm_error <- conditionMessage(fit)
  fit <- NULL
}
# the MLE under brms's draw names, through frm_sample()'s own map from
# the objective's coordinates to the draws
mle <- NULL
if (!is.null(fit)) {
  mle <- tryCatch({
    p <- fit$obj$env$last.par.best
    a <- array(p, c(1L, 1L, length(p)),
               dimnames = list(NULL, NULL, names(p)))
    m <- frmtmb.sample:::sample_finish(a, fit, list(laplace = FALSE))
    stats::setNames(as.numeric(m[1, ]), colnames(m))
  }, error = function(e) {
    rec$mle_error <<- conditionMessage(e)
    NULL
  })
}

# ---- brms code and data under default priors
bform <- eq$eq_eval(spec$bf_brms %||% spec$bf, "brms")
ft <- eq$eq_fam_text(spec, "brms")
family <- if (is.na(ft)) NULL else eq$eq_eval(ft, "brms")
prior <- eq$eq_eval(spec$prior_brms %||% spec$prior, "brms")
bargs <- c(list(bform, data = spec$data), spec$brms_args)
if (!is.null(family)) bargs$family <- family
if (!is.null(prior)) bargs$prior <- prior
t0 <- el()
code <- do.call(brms::make_stancode, bargs)
sdat <- suppressMessages(do.call(brms::make_standata, bargs))
rec$t_codegen <- el() - t0
mock <- function(sf) {
  suppressMessages(suppressWarnings(do.call(
    brms::brm, c(bargs, list(backend = "mock", mock_fit = sf,
                             rename = TRUE)))))
}
summ <- function(bf) {
  dr <- posterior::as_draws_array(bf)
  s <- posterior::summarise_draws(dr, "mean", "median", "sd", "mcse_mean",
                                  "rhat", "ess_bulk")
  as.data.frame(s)
}
diag_np <- function(np) {
  c(div = sum(np$Value[np$Parameter == "divergent__"]),
    maxtd = max(np$Value[np$Parameter == "treedepth__"]))
}

# ---- stanli
t0 <- el()
sm <- tryCatch(stanli_model(code = code, data = unclass(sdat)),
               error = function(e) e)
rec$t_stanli_setup <- el() - t0
if (inherits(sm, "error")) {
  rec$stanli_error <- conditionMessage(sm)
} else {
  t0 <- el()
  sfit <- tryCatch(sample_model(sm, chains = 4, warmup = 1000,
                                samples = 1000, seed = 1, refresh = 0,
                                parallel_chains = 1),
                   error = function(e) e)
  rec$t_stanli_sample <- el() - t0
  if (inherits(sfit, "error")) {
    rec$stanli_error <- conditionMessage(sfit)
  } else {
    t0 <- el()
    bf_s <- tryCatch(mock(as_stanfit(sfit)), error = function(e) e)
    rec$t_stanli_mock <- el() - t0
    if (inherits(bf_s, "error")) {
      rec$stanli_error <- paste("mock:", conditionMessage(bf_s))
    } else {
      rec$sum_stanli <- summ(bf_s)
      rec$np_stanli <- diag_np(brms::nuts_params(bf_s))
    }
  }
}
cat("stanli setup", rec$t_stanli_setup, "sample", rec$t_stanli_sample,
    "error", rec$stanli_error %||% "-", "\n")

# ---- rstan, compiled fresh
t0 <- el()
rm_ <- tryCatch(rstan::stan_model(model_code = code),
                error = function(e) e)
rec$t_rstan_compile <- el() - t0
if (inherits(rm_, "error")) {
  rec$rstan_error <- conditionMessage(rm_)
} else {
  t0 <- el()
  rfit <- tryCatch(suppressWarnings(rstan::sampling(
    rm_, data = sdat, chains = 4, warmup = 1000, iter = 2000, seed = 1,
    refresh = 0, cores = 1)), error = function(e) e)
  rec$t_rstan_sample <- el() - t0
  if (inherits(rfit, "error")) {
    rec$rstan_error <- conditionMessage(rfit)
  } else {
    bf_r <- tryCatch(mock(rfit), error = function(e) e)
    if (inherits(bf_r, "error")) {
      rec$rstan_error <- paste("mock:", conditionMessage(bf_r))
    } else {
      rec$sum_rstan <- summ(bf_r)
      rec$np_rstan <- diag_np(brms::nuts_params(bf_r))
    }
  }
}
cat("rstan compile", rec$t_rstan_compile, "sample", rec$t_rstan_sample,
    "error", rec$rstan_error %||% "-", "\n")

# ---- agreement
is_re <- function(v) grepl("^(r_|z_|zgp_|zs_|s_|Xme|err)", v) |
  grepl("^(lp__|lprior)$", v)
if (!is.null(rec$sum_stanli) && !is.null(rec$sum_rstan)) {
  a <- rec$sum_stanli
  b <- rec$sum_rstan[match(a$variable, rec$sum_rstan$variable), ]
  ok <- !is.na(b$variable)
  z <- (a$mean[ok] - b$mean[ok]) /
    sqrt(a$mcse_mean[ok]^2 + b$mcse_mean[ok]^2)
  names(z) <- a$variable[ok]
  z <- z[is.finite(z)]
  rec$z_sr <- z
  rec$sr_nvar <- length(z)
  rec$sr_maxz <- max(abs(z))
  rec$sr_frac3 <- mean(abs(z) > 3)
  rec$sr_crit <- stats::qnorm(1 - 0.025 / length(z))
  rec$sr_lp_z <- unname(z["lp__"])
  rec$sr_only_stanli <- setdiff(a$variable, rec$sum_rstan$variable)
  rec$sr_only_rstan <- setdiff(rec$sum_rstan$variable, a$variable)
}
ref <- rec$sum_stanli %||% rec$sum_rstan
if (!is.null(mle) && !is.null(ref)) {
  common <- intersect(names(mle), ref$variable)
  s <- ref[match(common, ref$variable), ]
  zm <- (mle[common] - s$median) / s$sd
  rec$z_mle <- zm
  hyp <- !is_re(common)
  rec$mle_nhyp <- sum(hyp)
  rec$mle_maxz_hyp <- if (any(hyp)) max(abs(zm[hyp])) else NA
  rec$mle_worst <- if (any(hyp)) common[hyp][which.max(abs(zm[hyp]))] else NA
  rec$mle_medz_re <- if (any(!hyp)) stats::median(abs(zm[!hyp])) else NA
  rec$mle_missing <- setdiff(names(mle), ref$variable)
}
for (k in c("stanli", "rstan")) {
  s <- rec[[paste0("sum_", k)]]
  if (!is.null(s)) {
    rec[[paste0(k, "_maxrhat")]] <- max(s$rhat, na.rm = TRUE)
    rec[[paste0(k, "_miness")]] <- min(s$ess_bulk, na.rm = TRUE)
  }
}
saveRDS(rec, file.path(outdir, paste0(id, ".rds")))
num <- function(x) if (is.null(x) || !length(x)) NA_real_ else
  round(as.numeric(x)[1], 3)
row <- data.frame(
  id = id, key = rec$key, n = rec$nobs,
  t_frm = num(rec$t_frm), t_codegen = num(rec$t_codegen),
  t_stanli_setup = num(rec$t_stanli_setup),
  t_stanli_sample = num(rec$t_stanli_sample),
  t_rstan_compile = num(rec$t_rstan_compile),
  t_rstan_sample = num(rec$t_rstan_sample),
  stanli_div = num(rec$np_stanli[["div"]]),
  rstan_div = num(rec$np_rstan[["div"]]),
  stanli_maxrhat = num(rec$stanli_maxrhat),
  rstan_maxrhat = num(rec$rstan_maxrhat),
  stanli_miness = num(rec$stanli_miness),
  rstan_miness = num(rec$rstan_miness),
  sr_nvar = num(rec$sr_nvar), sr_maxz = num(rec$sr_maxz),
  sr_crit = num(rec$sr_crit), sr_frac3 = num(rec$sr_frac3),
  sr_lp_z = num(rec$sr_lp_z),
  mle_nhyp = num(rec$mle_nhyp), mle_maxz_hyp = num(rec$mle_maxz_hyp),
  mle_worst = rec$mle_worst %||% NA, mle_medz_re = num(rec$mle_medz_re),
  frm_error = rec$frm_error %||% "", mle_error = rec$mle_error %||% "",
  stanli_error = substr(rec$stanli_error %||% "", 1, 300),
  rstan_error = substr(rec$rstan_error %||% "", 1, 300))
utils::write.table(row, file.path(outdir, paste0(id, ".tsv")), sep = "\t",
                   row.names = FALSE, quote = TRUE)
print(t(row))
