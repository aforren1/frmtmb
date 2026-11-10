# Reviewer scratch (review of lane translator, 2026-11-09): p048
# (t2() x rescor). A model mutant that scales t2() block 2's penalized
# part by 1.01 still read L1 TRUE (constants 3.3e-10, 5.2e-10, 5.1e-10;
# dev/review-translator-modmut.log). Which blocks does L1 see? Prints the
# three sds and the norm of each block's zs and s at the three points the
# survey reads, then runs map mutants that break ONE block each.
# Records nothing.
#
#   Rscript dev/review-translator-p048.R
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
fmt <- function(x) if (!length(x) || !is.numeric(x)) "-" else
  paste(signif(x, 3), collapse = ",")
fresh <- function() {
  e <- new.env(parent = asNamespace("frmtmb"))
  for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
    sys.source(f, envir = e)
  }
  e
}
smp <- read.delim("dev/equiv-data/sample.tsv")
id <- "p048"
row <- smp[smp$id == id, ]
key <- paste(row$feature_a, row$feature_b, sep = " x ")

# the translated blocks at the estimate and at the two stepped points
e <- fresh()
spec <- e$eq_specs[[key]]()
bform <- e$eq_eval(spec$bf, "brms")
fam <- brms::brmsfamily("gaussian")
prior <- e$eq_brms_prior(spec, bform, fam, spec$data)
args <- e$eq_brms_args(spec, spec$data, prior)
code <- do.call(brms::make_stancode, c(list(bform, family = fam), args))
sdat <- suppressMessages(do.call(brms::make_standata,
                                 c(list(bform, family = fam), args)))
rtab <- do.call(e$brms_ranef_table, c(list(bform, spec$data, fam, prior),
                                      spec$brms_args))
fits <- list(est = e$eq_frm_fit(spec))
for (k in 1:2) {
  fits[[paste0("step", k)]] <- e$eq_frm_fit(spec, ctrl = list(
    optimizer = e$eq_step_optimizer(c(1, -1)[k]), restarts = 0,
    check_se = "ignore"))
}
for (nm in names(fits)) {
  pars <- e$stan_pars_from_fit(fits[[nm]], sdat, code, rtab)
  sds <- as.numeric(pars[["sds_y1_1"]])
  zn <- vapply(1:3, function(j) {
    sqrt(sum(as.numeric(pars[[paste0("zs_y1_1_", j)]])^2))
  }, 0)
  cat(sprintf("%-6s sds_y1_1 %s | |zs_j| %s | |s_j| = sds*|zs| %s\n", nm,
              fmt(sds), fmt(zn), fmt(sds * zn)))
}

run <- function(label, mutate = NULL, mut_code = NULL) {
  e <- fresh()
  if (is.function(mutate)) mutate(e)
  if (is.function(mut_code)) {
    orig <- e$eq_l1
    e$eq_l1 <- function(..., code_fix = NULL) {
      orig(..., code_fix = mut_code)
    }
  }
  r <- e$eq_measure(key, row)
  cat(sprintf("%-44s verdict %-12s L1 ok %-5s const %s grad %s err %s\n",
              label, r$verdict %||% "", format(r$l1_ok %||% NA),
              fmt(c(r$l1_const, r$l1_extra_const)), fmt(r$l1_grad),
              substr(r$l1_error %||% "", 1, 100)))
}
for (j in 1:3) {
  run(sprintf("map: zs_y1_1_%d alone reversed", j), mutate = function(e) {
    v <- e$eq_rules$zs_prefixed$value
    e$eq_rules$zs_prefixed$value <- function(nm, ...) {
      out <- v(nm, ...)
      if (identical(nm, paste0("zs_y1_1_", j))) {
        out$value <- array(rev(out$value), length(out$value))
      }
      out
    }
  })
}
for (j in c(1, 3)) {
  from <- sprintf("Zs_y1_1_%d * s_y1_1_%d", j, j)
  to <- sprintf("Zs_y1_1_%d * (1.01 * s_y1_1_%d)", j, j)
  run(sprintf("model: t2 block %d scaled 1.01", j), mut_code = local({
    f <- from
    t <- to
    function(code) {
      new <- sub(f, t, code, fixed = TRUE)
      if (identical(new, code)) stop("pattern not found")
      new
    }
  }))
}
