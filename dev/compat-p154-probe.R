# p154 (cens() x cosy()) read "refusal correct" while its row was
# refused. The live row is now conditional: refused on gaussian and
# student, working on a family with no residual (brms's latent residual).
# The sampled spec is gaussian, so it does not meet the condition, and
# scored against the live row it reads "compat row wrong". This probe
# runs the flow on the same pair spelled on poisson, where the condition
# holds, without editing dev/equiv-specs.R.
#   Rscript dev/compat-p154-probe.R <lib;lib>
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
eq_env$eq_specs[[key]] <- function() {
  d <- eq_env$eq_dat("poisson", seed = 154)
  d$c <- ifelse(d$y > 4, "right", "none")
  d$y <- pmin(d$y, 4)
  list(bf = "bf(y | cens(c) ~ x + cosy(gr = g))", family = "poisson()",
       data = d)
}
smp <- read.delim("dev/equiv-data/sample.tsv")
row <- smp[smp$id == "p154", ]
row$status <- "conditional"
rec <- eq_env$eq_measure(key, row)
cat("censored rows:", sum(eq_env$eq_specs[[key]]()$data$c == "right"),
    "of 200\n")
cat("verdict:", rec$verdict, "| level:", paste(rec$level, collapse = "+"),
    "| reason:", rec$reason, "\n")
cat("frm_error:", rec$frm_error %||% "<none>", "\n")
cat("brms_error:", rec$brms_error %||% "<none>", "\n")
cat("l1 const:", format(rec$l1_const, digits = 6), "grad:",
    format(rec$l1_grad, digits = 3), "points:", rec$l1_points, "ok:",
    rec$l1_ok, "extra:", format(rec$l1_extra_const, digits = 6), "\n")
