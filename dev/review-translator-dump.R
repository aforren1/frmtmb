# Reviewer scratch: write brms 2.23.0's program for each reviewed pair to
# dev/review-translator-out/code-<id>.stan (make_stancode, no compile).
.libPaths(c("C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); requireNamespace("brms")})
attach(list(me = brms::me), name = "equiv_brms_me", warn.conflicts = FALSE)
eq_env <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = eq_env)
}
dir.create("dev/review-translator-out", showWarnings = FALSE)
smp <- read.delim("dev/equiv-data/sample.tsv")
for (id in c("p042", "p048", "p131", "p021", "p044", "p136")) {
  row <- smp[smp$id == id, ]
  key <- paste(row$feature_a, row$feature_b, sep = " x ")
  spec <- eq_env$eq_specs[[key]]()
  bform <- eq_env$eq_eval(spec$bf, "brms")
  fam <- eq_env$eq_fam_text(spec, "brms")
  fam <- if (is.na(fam)) brms::brmsfamily("gaussian") else
    eq_env$eq_eval(fam, "brms")
  pr <- eq_env$eq_brms_prior(spec, bform, fam, spec$data)
  code <- do.call(brms::make_stancode,
                  c(list(bform, family = fam),
                    eq_env$eq_brms_args(spec, spec$data, pr)))
  writeLines(code, sprintf("dev/review-translator-out/code-%s.stan", id))
  cat(id, key, "written\n")
}
