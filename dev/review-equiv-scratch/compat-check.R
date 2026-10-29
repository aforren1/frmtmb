# Reviewer: spot-check compat-row proposals by running both packages on
# the lane's spec, and print the current frm_compat() row.
.libPaths(c("C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/review-equiv-cache",
                                             winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
attach(list(me = brms::me), name = "equiv_brms_me", warn.conflicts = FALSE)
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
keys <- c("cens() x disc", "cens() x equidistant", "vint() x sum_to_zero",
          "trials() x cs_pred()", "nbinom1 x nl", "zero_inflated_beta x nl",
          "cs_pred() x nl", "us x importance")
for (key in keys) {
  ft <- strsplit(key, " x ", fixed = TRUE)[[1]]
  st <- as.data.frame(frm_compat(ft[1], ft[2]))
  spec <- E$eq_specs[[key]]()
  if (length(spec$data2)) spec$brms_args$data2 <- spec$data2
  r <- E$eq_try(E$eq_frm_fit(spec))
  b <- if (!is.null(spec$nobrms)) list(error = "no brms analogue") else
    E$eq_try({
      bform <- E$eq_eval(spec$bf, "brms")
      fam <- E$eq_eval(spec$family %||% "gaussian()", "brms")
      pr <- E$eq_brms_prior(spec, bform, fam, spec$data)
      code <- do.call(brms::make_stancode,
                      c(list(bform, family = fam),
                        E$eq_brms_args(spec, spec$data, pr)))
      E$brms_stan_model(code)
      TRUE
    })
  cat(sprintf("== %s | compat: %s\n   frmtmb: %s\n   brms:   %s\n", key,
              st$status[1],
              if (is.null(r$error)) paste("fits, logLik",
                round(as.numeric(logLik(r$value)), 3)) else
                substr(r$error, 1, 160),
              if (is.null(b$error)) "builds and compiles" else
                substr(b$error, 1, 160)))
}
