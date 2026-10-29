# Repro for the pairs where frmtmb refuses a model brms 2.23.0 builds
# AND compiles (the survey confirms brms by one cached compile before it
# calls a refusal a divergence): |ID| x mmc() (p070), gr_cov x mmc()
# (p087), s() x prior with class "sds" (p043), von_mises x mi() (p169).
# Each prints frmtmb's error and whether brms's program compiles.
.libPaths(c("C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/review-equiv-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
for (key in c("|ID| x mmc()", "gr_cov x mmc()", "s() x prior",
              "von_mises x mi()")) {
  spec <- E$eq_specs[[key]]()
  if (length(spec$data2)) spec$brms_args$data2 <- spec$data2
  cat("==", key, "\n")
  r <- E$eq_try(E$eq_frm_fit(spec))
  cat("frmtmb:", if (is.null(r$error)) "fits" else
    substr(r$error, 1, 300), "\n")
  b <- E$eq_try({
    bform <- E$eq_eval(spec$bf, "brms")
    fam <- E$eq_eval(spec$family %||% "gaussian()", "brms")
    pr <- E$eq_brms_prior(spec, bform, fam, spec$data)
    code <- do.call(brms::make_stancode, c(list(bform, family = fam),
                                           E$eq_brms_args(spec, spec$data,
                                                          pr)))
    E$brms_stan_model(code)
    TRUE
  })
  cat("brms 2.23.0:", if (is.null(b$error)) "builds and compiles" else
    substr(b$error, 1, 300), "\n")
}
