# Question 2: could stanli stand in for rstan in the equivalence survey
# (dev/equiv-measure.R) and the brms test tiers?
#
#   Rscript dev/samplers-q2.R <pNNN> [outdir]
#
# The survey reads brms's program in two ways:
#   L1  rstan::log_prob() / grad_log_prob() WITHOUT the transform
#       Jacobian (adjust_transform = FALSE) at frm()'s estimate, and the
#       Jacobian of the prior blocks; on a stanfit from
#       rstan::sampling(chains = 0) of a compiled (cached) program.
#   L4  rstan::sampling(), 2 chains x 2000, flat priors.
# stanli's density always includes the Jacobian (its README). The shim
# measured here: a second stanli program, brms's code with its model
# block emptied, whose log density is the Jacobian alone; the density
# without the Jacobian is the difference. It is compared with rstan's
# adjust_transform = FALSE at frm()'s estimate and at 5 points around it
# (seed 1). The times of each Stan step are measured on both sides.
args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
outdir <- if (length(args) >= 2) args[2] else "dev/samplers-q2-out"
source("dev/samplers-designs.R")
samplers_lib()
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache",
                                             winslash = "/"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
eq <- samplers_eq_env()
suppressPackageStartupMessages({
  library(stanli)
  library(posterior)
})
el <- function() proc.time()[["elapsed"]]
rec <- list(id = id)
spec <- samplers_design(id, eq)
rec$key <- spec$key
cat("pair", id, ":", spec$key, "\n")
fit <- eq$eq_frm_fit(spec)
bform <- eq$eq_eval(spec$bf_brms %||% spec$bf, "brms")
ft <- eq$eq_fam_text(spec, "brms")
family <- if (is.na(ft)) NULL else eq$eq_eval(ft, "brms")
prior <- eq$eq_brms_prior(spec, bform, family, spec$data)
data <- eq$brms_mi_observed(fit, spec$data)
bargs <- c(list(bform, family = family), eq$eq_brms_args(spec, data, prior))
code <- do.call(brms::make_stancode, bargs)
sdat <- suppressMessages(do.call(brms::make_standata, bargs))
rtab <- do.call(eq$brms_ranef_table,
                c(list(bform, data, family, prior), spec$brms_args))

# ---- rstan: warm cache, cold compile, chains = 0 stanfit
key <- eq$brms_stan_cache_key(code)
cache <- file.path(Sys.getenv("FRMTMB_STAN_CACHE"), paste0(key, ".rds"))
rec$cache_hit <- file.exists(cache)
if (rec$cache_hit) {
  t0 <- el()
  mod <- readRDS(cache)
  rec$t_rstan_warm <- el() - t0
}
t0 <- el()
mod_cold <- rstan::stan_model(model_code = code)
rec$t_rstan_compile <- el() - t0
if (!rec$cache_hit) mod <- mod_cold
t0 <- el()
sf <- suppressMessages(rstan::sampling(mod, data = sdat, chains = 0))
rec$t_rstan_chains0 <- el() - t0
pars <- eq$stan_pars_from_fit(fit, sdat, code, rtab)
u0 <- rstan::unconstrain_pars(sf, pars)
set.seed(1)
U <- rbind(u0, t(replicate(5, u0 + rnorm(length(u0), 0, 0.3))))
t0 <- el()
r_noj <- apply(U, 1, function(u) rstan::log_prob(sf, u,
                                                 adjust_transform = FALSE))
r_j <- apply(U, 1, function(u) rstan::log_prob(sf, u,
                                               adjust_transform = TRUE))
r_g <- t(apply(U, 1, function(u) rstan::grad_log_prob(
  sf, u, adjust_transform = FALSE)))
rec$t_rstan_lp <- el() - t0

# ---- stanli: the program and its Jacobian-only twin
jac_code <- function(code) {
  lines <- strsplit(code, "\n", fixed = TRUE)[[1]]
  m0 <- grep("^model \\{", lines)
  gq <- grep("^generated quantities \\{", lines)
  end <- if (length(gq)) gq - 1L else length(lines)
  # the model block's closing brace is the last "}" before `end`
  close <- max(which(lines[seq_len(end)] == "}"))
  paste(c(lines[seq_len(m0)], lines[(close):length(lines)]),
        collapse = "\n")
}
t0 <- el()
sm <- stanli_model(code = code, data = unclass(sdat))
rec$t_stanli_setup <- el() - t0
t0 <- el()
smj <- stanli_model(code = jac_code(code), data = unclass(sdat))
rec$t_stanli_setup_jac <- el() - t0
# log_prob_grad()'s value is read by position or attribute, whichever
# this stanli build returns
lpg <- function(m, u) {
  r <- log_prob_grad(m, u)
  if (is.list(r)) list(lp = as.numeric(r[[1]]), g = as.numeric(r[[2]])) else
    list(lp = as.numeric(r[1]), g = as.numeric(attr(r, "gradient")))
}
rec$lpg_shape <- utils::capture.output(str(log_prob_grad(sm, U[1, ])))
t0 <- el()
s_j <- apply(U, 1, function(u) lpg(sm, u)$lp)
s_jac <- apply(U, 1, function(u) lpg(smj, u)$lp)
s_g <- t(apply(U, 1, function(u) lpg(sm, u)$g - lpg(smj, u)$g))
rec$t_stanli_lp <- el() - t0
s_noj <- s_j - s_jac
rec$lp <- data.frame(rstan_noj = r_noj, stanli_noj = s_noj,
                     rstan_j = r_j, stanli_j = s_j)
rec$d_noj <- max(abs(s_noj - r_noj) / pmax(1, abs(r_noj)))
rec$d_j <- max(abs(s_j - r_j) / pmax(1, abs(r_j)))
rec$d_grad <- max(abs(s_g - r_g) / pmax(1, abs(r_g)))
# the unconstrained parameterization is the same: stanli's unconstrain()
# of rstan's constrained point returns rstan's vector
uq <- tryCatch(unconstrain(sm, rstan::constrain_pars(sf, u0),
                           ignore_unknown = TRUE),
               error = function(e) e)
rec$d_unconstrain <- if (inherits(uq, "error")) conditionMessage(uq) else
  max(abs(as.numeric(uq) - u0))
cat("L1 shim: rel diff no-Jacobian", rec$d_noj, "with", rec$d_j,
    "gradient", rec$d_grad, "unconstrain", rec$d_unconstrain, "\n")

# ---- L4: 2 chains x 2000 with flat priors, both samplers
t0 <- el()
rfit <- suppressWarnings(rstan::sampling(mod, data = sdat, chains = 2,
                                         iter = 2000, seed = 1, refresh = 0,
                                         cores = 1))
rec$t_rstan_l4 <- el() - t0
t0 <- el()
sfit <- sample_model(sm, chains = 2, warmup = 1000, samples = 1000,
                     seed = 1, refresh = 0, parallel_chains = 1)
rec$t_stanli_l4 <- el() - t0
sa <- as.data.frame(summarise_draws(as_draws_array(rfit), "mean",
                                    "mcse_mean", "rhat"))
sb <- as.data.frame(summarise_draws(as_draws_array(as_stanfit(sfit)),
                                    "mean", "mcse_mean", "rhat"))
sb <- sb[match(sa$variable, sb$variable), ]
ok <- !is.na(sb$variable)
z <- (sa$mean[ok] - sb$mean[ok]) / sqrt(sa$mcse_mean[ok]^2 +
                                          sb$mcse_mean[ok]^2)
z <- z[is.finite(z)]
rec$l4_maxz <- max(abs(z))
rec$l4_nvar <- length(z)
rec$l4_rhat <- c(rstan = max(sa$rhat, na.rm = TRUE),
                 stanli = max(sb$rhat, na.rm = TRUE))
# ---- the rstan-shaped calls the test tiers make, on as_stanmodel()
asm <- as_stanmodel(code)
z0 <- tryCatch(rstan::sampling(asm, data = sdat, chains = 0),
               error = function(e) e)
rec$asm_chains0 <- if (inherits(z0, "error")) conditionMessage(z0) else
  class(z0)[1]
if (!inherits(z0, "error")) {
  e <- tryCatch(rstan::log_prob(z0, u0, adjust_transform = FALSE),
                error = function(e) e)
  rec$asm_noj <- if (inherits(e, "error")) conditionMessage(e) else e
}
# the method tier's Fixed_param brmsfit (helper-brms-methods.R)
scaffold <- suppressMessages(do.call(brms::brm, c(bargs, list(empty = TRUE))))
excl <- get("exclude_pars", asNamespace("brms"))(scaffold)
fp <- function(model) {
  sc <- scaffold
  sc$fit <- suppressMessages(rstan::sampling(
    model, data = sdat, chains = 1, iter = 10, warmup = 0,
    algorithm = "Fixed_param", init = list(pars), refresh = 0, seed = 1,
    pars = excl, include = FALSE))
  brms::rename_pars(sc)
}
t0 <- el()
fr <- tryCatch(fp(mod), error = function(e) e)
rec$t_fp_rstan <- el() - t0
t0 <- el()
fs <- tryCatch(fp(asm), error = function(e) e)
rec$t_fp_stanli <- el() - t0
if (inherits(fs, "error")) {
  rec$fp_error <- conditionMessage(fs)
} else if (!inherits(fr, "error")) {
  ea <- brms::posterior_epred(fr)
  eb <- brms::posterior_epred(fs)
  rec$fp_epred_diff <- max(abs(ea - eb) / pmax(1, abs(ea)))
  rec$fp_vars_same <- identical(sort(posterior::variables(fr)),
                                sort(posterior::variables(fs)))
}
cat("as_stanmodel chains = 0:", rec$asm_chains0, "| no-Jacobian log_prob:",
    format(rec$asm_noj), "| Fixed_param:", rec$fp_error %||% "ok",
    rec$fp_epred_diff %||% "", rec$fp_vars_same %||% "", "\n")

saveRDS(rec, file.path(outdir, paste0(id, ".rds")))
num <- function(x) if (is.null(x) || !length(x) || is.character(x))
  NA_real_ else signif(as.numeric(x)[1], 4)
row <- data.frame(
  id = id, key = rec$key, npar_u = length(u0), cache_hit = rec$cache_hit,
  t_rstan_warm = num(rec$t_rstan_warm),
  t_rstan_compile = num(rec$t_rstan_compile),
  t_rstan_chains0 = num(rec$t_rstan_chains0),
  t_rstan_lp = num(rec$t_rstan_lp),
  t_stanli_setup = num(rec$t_stanli_setup),
  t_stanli_setup_jac = num(rec$t_stanli_setup_jac),
  t_stanli_lp = num(rec$t_stanli_lp),
  d_noj = num(rec$d_noj), d_j = num(rec$d_j), d_grad = num(rec$d_grad),
  d_unconstrain = num(rec$d_unconstrain),
  t_rstan_l4 = num(rec$t_rstan_l4), t_stanli_l4 = num(rec$t_stanli_l4),
  l4_nvar = rec$l4_nvar, l4_maxz = num(rec$l4_maxz),
  l4_crit = num(stats::qnorm(1 - 0.025 / rec$l4_nvar)),
  l4_rhat_rstan = num(rec$l4_rhat[["rstan"]]),
  l4_rhat_stanli = num(rec$l4_rhat[["stanli"]]),
  fp_epred_diff = num(rec$fp_epred_diff),
  fp_vars_same = rec$fp_vars_same %||% NA,
  t_fp_rstan = num(rec$t_fp_rstan), t_fp_stanli = num(rec$t_fp_stanli),
  asm_chains0 = substr(rec$asm_chains0, 1, 120),
  asm_noj = substr(format(rec$asm_noj), 1, 160),
  fp_error = substr(rec$fp_error %||% "", 1, 200))
utils::write.table(row, file.path(outdir, paste0(id, ".tsv")), sep = "\t",
                   row.names = FALSE, quote = TRUE)
print(t(row))
