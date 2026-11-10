# Reviewer scratch (review of lane translator, 2026-11-09): MODEL mutants.
# The lane's mutants break the map; these break the model. brms's program
# is edited by a fixed string substitution (the map is left alone), so
# the two packages' models differ in one place. Every mutant must turn
# its pair's L1 (or, for p136, the patched-program L1 of the decision
# path) from pass to fail. Records nothing.
#
#   Rscript dev/review-translator-modmut.R
lib <- "C:/Users/adf44/source/r/rellib-r19"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache",
                                             winslash = "/"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
  requireNamespace("brms")
  requireNamespace("rstan")
})
attach(list(me = brms::me), name = "equiv_brms_me", warn.conflicts = FALSE)
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
`%||%` <- function(a, b) if (is.null(a)) b else a
fmt <- function(x) if (!length(x) || !is.numeric(x)) "-" else paste(signif(x, 3), collapse = ",")
fresh <- function() {
  e <- new.env(parent = asNamespace("frmtmb"))
  for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
    sys.source(f, envir = e)
  }
  e
}
smp <- read.delim("dev/equiv-data/sample.tsv")
sub1 <- function(from, to) {
  force(from)
  force(to)
  function(code) {
    new <- sub(from, to, code, fixed = TRUE)
    if (identical(new, code)) stop("mutant pattern not found: ", from)
    new
  }
}
run <- function(label, id, mut = NULL) {
  e <- fresh()
  if (is.function(mut)) {
    orig <- e$eq_l1
    e$eq_l1 <- function(..., code_fix = NULL) {
      cf <- code_fix
      orig(..., code_fix = function(code) {
        mut(if (is.function(cf)) cf(code) else code)
      })
    }
  }
  row <- smp[smp$id == id, ]
  key <- paste(row$feature_a, row$feature_b, sep = " x ")
  t0 <- proc.time()[["elapsed"]]
  r <- e$eq_measure(key, row)
  cat(sprintf(paste0("%-40s %s  verdict %-22s level %-6s L1 ok %-5s ",
                     "const %s  dec %s dec_ok %s  err %s  (%.0f s)\n"),
              label, id, r$verdict %||% "", r$level %||% "",
              format(r$l1_ok %||% NA),
              fmt(c(r$l1_const, r$l1_extra_const)),
              fmt(c(r$l1_dec_const, r$l1_dec_grad)),
              format(r$l1_dec_ok %||% NA),
              substr(paste(r$l1_error %||% "", r$l1_dec_error %||% ""), 1,
                     120),
              proc.time()[["elapsed"]] - t0))
  invisible(r)
}

# smooth in an nl body: the smooth's penalized part scaled by 1.01 (a
# different penalty scale), and the unpenalized columns dropped
run("p042 smooth scaled 1.01", "p042",
    sub1("Zs_b1_1_1 * s_b1_1_1", "Zs_b1_1_1 * (1.01 * s_b1_1_1)"))
run("p042 unpenalized part dropped", "p042",
    sub1("+ Xs_b1 * bs_b1 ", "+ 0 * Xs_b1 * bs_b1 "))
# t2() under mvbf: the second penalty block scaled by 1.01
run("p048 t2 block 2 scaled 1.01", "p048",
    sub1("Zs_y1_1_2 * s_y1_1_2", "Zs_y1_1_2 * (1.01 * s_y1_1_2)"))
# sratio mixture with cs(): the cs() sign of component 2, and its disc
run("p131 cs sign of component 2", "p131",
    sub1("Intercept_mu2 - transpose(mucs_mu2[n])",
         "Intercept_mu2 + transpose(mucs_mu2[n])"))
run("p131 disc2 = 1.1", "p131",
    sub1("real disc2 = 1;", "real disc2 = 1.1;"))
# gev mixture: the shape enters the density at 1.02 times its value
run("p021 gev t = 1 + 1.02 xi x", "p021",
    sub1("real t = 1 + xi * x;", "real t = 1 + 1.02 * xi * x;"))
# gp in an nl body: a length scale 1.01 times brms's, in BOTH programs
# (the plain and the nugget-patched one); must not score "by decision"
run("p136 lscale 1.01 (both programs)", "p136",
    sub1("cov = gp_exp_quad_cov(x, sdgp, lscale[1]);",
         "cov = gp_exp_quad_cov(x, sdgp, 1.01 * lscale[1]);"))
# and a nugget 100 times brms's: the patch no longer matches its line
run("p136 brms nugget 1e-10", "p136",
    function(code) {
      new <- gsub("1e-12", "1e-10", code, fixed = TRUE)
      if (identical(new, code)) stop("no nugget")
      new
    })
