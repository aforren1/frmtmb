# Probe one survey spec: what brms declares, what frmtmb's fit holds,
# and where the translator stops. A diagnostic for lane translator; it
# records nothing.
#
#   Rscript dev/translator-probe.R "<feature_a> x <feature_b>"
#
# Run from the worktree root with FRMTMB_EQUIV_LIB set as for
# dev/equiv-measure.R.
args <- commandArgs(trailingOnly = TRUE)
key <- args[1]
lib <- Sys.getenv("FRMTMB_EQUIV_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE"))) {
  Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache",
                                               winslash = "/"))
}
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
  requireNamespace("brms")
  requireNamespace("rstan")
})
attach(list(me = brms::me), name = "equiv_brms_me", warn.conflicts = FALSE)
eq_env <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = eq_env)
}
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
probe <- local({
  spec <- eq_specs[[key]]()
  if (length(spec$data2)) spec$brms_args$data2 <- spec$data2
  fit <- eq_frm_fit(spec)
  bform <- eq_eval(spec$bf_brms %||% spec$bf, "brms")
  ft <- eq_fam_text(spec, "brms")
  family <- if (is.na(ft)) NULL else eq_eval(ft, "brms")
  prior <- eq_brms_prior(spec, bform, family, spec$data)
  args <- eq_brms_args(spec, spec$data, prior)
  code <- do.call(brms::make_stancode, c(list(bform, family = family), args))
  sdat <- suppressMessages(do.call(brms::make_standata,
                                   c(list(bform, family = family), args)))
  rtab <- tryCatch(do.call(brms_ranef_table,
                           c(list(bform, spec$data, family, prior),
                             spec$brms_args)), error = function(e) NULL)
  list(spec = spec, fit = fit, code = code, sdat = sdat, rtab = rtab,
       bform = bform, family = family, prior = prior)
}, envir = eq_env)
assign("probe", probe, envir = globalenv())
fit <- probe$fit
cat("\n== Stan parameters\n")
print(eq_env$brms_stan_par_names(probe$code))
cat("\n== standata names\n")
print(names(probe$sdat))
cat("\n== frmtmb estimates\n")
str(fit$estimates, give.attr = FALSE)
cat("\n== fixef_by_dpar names\n")
print(lapply(frmtmb:::fixef_by_dpar(fit), names))
cat("\n== linpreds\n")
print(names(fit$frame[["linpreds"]]))
cat("\n== re_blocks\n")
for (b in fit$frame[["re_blocks"]]) {
  cat(sprintf("  cs=%s dpar=%s label=%s lp_key=%s dim=%s\n",
              b[["covstruct"]], b[["dpar"]] %||% "<NULL>",
              b[["term_label"]] %||% "", b$components[[1]]$lp_key %||% "",
              b[["dim"]] %||% ""))
}
cat("\n== family links\n")
print(names(family(fit)$links))
cat("\n== translator\n")
r <- tryCatch(withCallingHandlers(
  eq_env$stan_pars_from_fit(fit, probe$sdat, probe$code, probe$rtab),
  error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    cs <- sys.calls()
    for (k in seq_along(cs)) cat(k, ":", deparse1(cs[[k]])[1], "\n")
  }), error = function(e) NULL)
if (!is.null(r)) str(r)
if (!is.null(probe$sdat[["Xs"]])) {
  cat("\n== log|det| of mu's bs basis map:",
      format(eq_env$eq_bs_logdet(fit, probe$sdat), digits = 17), "\n")
}
