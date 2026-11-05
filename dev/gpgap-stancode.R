# Print brms 2.23.0's Stan code and standata for the exact gp() pairs,
# the gp() pieces only. No compile.
#   Rscript dev/gpgap-stancode.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
keys <- c("lognormal x gp", "gp x nl")
local(envir = E, {
  for (key in keys) {
    spec <- eq_specs[[key]]()
    bform <- eq_eval(spec$bf, "brms")
    ft <- eq_fam_text(spec, "brms")
    family <- if (is.na(ft)) NULL else eq_eval(ft, "brms")
    prior <- eq_brms_prior(spec, bform, family, spec$data)
    args <- eq_brms_args(spec, spec$data, prior)
    code <- do.call(brms::make_stancode, c(list(bform, family = family),
                                            args))
    sdat <- suppressMessages(do.call(brms::make_standata,
                                     c(list(bform, family = family), args)))
    cat("=====", key, "\n")
    cat(code, "\n")
    g <- grep("gp", names(sdat), value = TRUE)
    for (nm in g) {
      cat("--", nm, ": "); str(sdat[[nm]])
    }
  }
})
