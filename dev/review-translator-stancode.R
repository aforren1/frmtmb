# Reviewer scratch: brms 2.23.0's program for p021 (gev mixture) and the
# gp nugget line of p136, read with make_stancode (no compile).
.libPaths(c("C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); requireNamespace("brms")})
cat("brms", as.character(packageVersion("brms")), "\n")
eq_env <- new.env(parent = asNamespace("frmtmb"))
sys.source("tests/testthat/helper-brms.R", envir = eq_env)
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = eq_env)
}
for (key in c("gen_extreme_value x mixture", "gp x nl")) {
  spec <- eq_env$eq_specs[[key]]()
  bform <- eval(parse(text = spec$bf), envir = asNamespace("brms"))
  fam <- if (is.null(spec$family)) stats::gaussian() else
    eval(parse(text = spec$family), envir = asNamespace("brms"))
  pr <- eq_env$brms_flat_prior(bform, data = spec$data, family = fam)
  code <- brms::make_stancode(bform, data = spec$data, family = fam,
                              prior = pr)
  cat("=====", key, "\n")
  ln <- strsplit(code, "\n")[[1]]
  hit <- grep("xi|1e-12|cov\\[|gp_exp|sdgp|lprior|target", ln)
  cat(paste(hit, ln[hit]), sep = "\n")
}
