# Question 2 follow-up: the method tier's Fixed_param brmsfit
# (tests/testthat/helper-brms-methods.R) on stanli's as_stanmodel().
# dev/samplers-q2.R found that stanli refuses brms's exclusion list
# because it names variables the program lacks, which rstan ignores.
# Cut to the program's own variables, the list is accepted, but stanli
# 0.25.0 then drops EVERY variable under include = FALSE (only lp__
# is left; dev/samplers-q2fp-debug.R). So stanli is given the complement,
# the variables to keep, with include = TRUE (two lines a swap would
# add), and the two brmsfits are compared.
#
#   Rscript dev/samplers-q2fp.R <pNNN> [<pNNN> ...]
ids <- commandArgs(trailingOnly = TRUE)
source("dev/samplers-designs.R")
samplers_lib()
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache",
                                             winslash = "/"))
eq <- samplers_eq_env()
suppressPackageStartupMessages(library(stanli))
el <- function() proc.time()[["elapsed"]]
out <- list()
for (id in ids) {
  spec <- samplers_design(id, eq)
  fit <- eq$eq_frm_fit(spec)
  bform <- eq$eq_eval(spec$bf_brms %||% spec$bf, "brms")
  ft <- eq$eq_fam_text(spec, "brms")
  family <- if (is.na(ft)) NULL else eq$eq_eval(ft, "brms")
  prior <- eq$eq_brms_prior(spec, bform, family, spec$data)
  bargs <- c(list(bform, family = family),
             eq$eq_brms_args(spec, spec$data, prior))
  code <- do.call(brms::make_stancode, bargs)
  sdat <- suppressMessages(do.call(brms::make_standata, bargs))
  scaffold <- suppressMessages(do.call(brms::brm,
                                       c(bargs, list(empty = TRUE))))
  pars <- eq$stan_pars_from_fit(fit, sdat, code, scaffold$ranef)
  mod <- eq$brms_stan_model(code)
  sf0 <- suppressMessages(rstan::sampling(mod, data = sdat, chains = 0))
  excl <- get("exclude_pars", asNamespace("brms"))(scaffold)
  excl_cut <- intersect(excl, sf0@model_pars)
  keep <- setdiff(sf0@model_pars, c(excl, "lp__"))
  fp <- function(model, ex, include = FALSE) {
    sc <- scaffold
    sc$fit <- suppressMessages(rstan::sampling(
      model, data = sdat, chains = 1, iter = 10, warmup = 0,
      algorithm = "Fixed_param", init = list(pars), refresh = 0, seed = 1,
      pars = ex, include = include))
    brms::rename_pars(sc)
  }
  t0 <- el()
  fr <- fp(mod, excl)
  t_r <- el() - t0
  t0 <- el()
  asm <- as_stanmodel(code)
  fs <- tryCatch(fp(asm, keep, include = TRUE), error = function(e) e)
  t_s <- el() - t0
  row <- data.frame(id = id, key = spec$key, n_excl = length(excl),
                    n_excl_cut = length(excl_cut),
                    t_fp_rstan = round(t_r, 3), t_fp_stanli = round(t_s, 3))
  if (inherits(fs, "error")) {
    row$error <- conditionMessage(fs)
  } else {
    ea <- brms::posterior_epred(fr)
    eb <- brms::posterior_epred(fs)
    row$epred_rel_diff <- signif(max(abs(ea - eb) / pmax(1, abs(ea))), 3)
    va <- posterior::variables(fr)
    vb <- posterior::variables(fs)
    row$vars_same <- identical(sort(va), sort(vb))
    row$fixef_rel_diff <- signif(max(abs(brms::fixef(fr) -
                                           brms::fixef(fs))), 3)
    row$error <- ""
  }
  print(row)
  out[[id]] <- row
}
dir.create("dev/samplers-q2-out", showWarnings = FALSE)
utils::write.table(do.call(rbind, out), "dev/samplers-q2-out/fp.tsv.txt",
                   sep = "\t", row.names = FALSE, quote = TRUE)
