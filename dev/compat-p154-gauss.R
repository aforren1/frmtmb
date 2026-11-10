# p154 (cens() x cosy()) in its old gaussian spelling, scored by the
# edited flow against the live row (conditional): both packages refuse,
# which is the refused half of the condition, so the verdict must be
# "refusal correct" with no proposal (review 2026-11-10, m2).
#   Rscript dev/compat-p154-gauss.R <lib;lib>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(args[1], ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache",
                                             winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
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
key <- "cens() x cosy()"
# the spec as it stood at e0cce906
eq_env$eq_specs[[key]] <- function() {
  d <- eq_env$eq_dat(seed = 154)
  d$c <- ifelse(d$y > 1.5, "right", "none")
  d$y <- pmin(d$y, 1.5)
  list(bf = "bf(y | cens(c) ~ x + cosy(gr = g))", data = d)
}
smp <- read.delim("dev/equiv-data/sample.tsv")
row <- smp[smp$id == "p154", ]
for (st in c("conditional", "refused", "works")) {
  row$status <- st
  rec <- eq_env$eq_measure(key, row)
  cat("live", st, "-> verdict:", rec$verdict, "| reason:", rec$reason,
      "| proposal:", if (is.null(rec$proposal)) "<none>" else rec$proposal,
      "\n")
}
cat("frm_error:", substr(rec$frm_error %||% "<none>", 1, 90), "\n")
cat("brms_error:", rec$brms_error %||% "<none>", "\n")
