# Reviewer scratch (review of lane translator, 2026-11-09): rule branches
# the lane's pairs never exercise.
#   1. ord_mix_intercept's centering shift: p131's components have only
#      cs() terms, so X_mu<k> is empty and the shift is 0. A variant with
#      a plain x in each component reads the shift; a mutant flips its
#      sign.
#   2. bcs_prefixed's byrow layout: p131 has one cs() column, where byrow
#      cannot matter. A variant with cs(x) + cs(z); a mutant reads
#      bycol.
#   3. eq_bs_logdet(): 0 on p044. A REML variant with a 2-D smooth;
#      a mutant drops the log-determinant.
# Records nothing.
#
#   Rscript dev/review-translator-branches.R
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
run <- function(label, id, mutate = NULL, spec_edit = NULL) {
  e <- fresh()
  if (is.function(mutate)) mutate(e)
  row <- smp[smp$id == id, ]
  key <- paste(row$feature_a, row$feature_b, sep = " x ")
  if (is.function(spec_edit)) {
    mk <- e$eq_specs[[key]]
    e$eq_specs[[key]] <- function() spec_edit(mk())
  }
  t0 <- proc.time()[["elapsed"]]
  r <- e$eq_measure(key, row)
  cat(sprintf(paste0("%-40s %s  verdict %-14s level %-6s L1 ok %-5s ",
                     "const %s  marg %s  err %s  (%.0f s)\n"),
              label, id, r$verdict %||% "", r$level %||% "",
              format(r$l1_ok %||% NA),
              fmt(c(r$l1_const, r$l1_extra_const)),
              fmt(r$l1_marg),
              substr(paste(r$l1_error %||% "", r$frm_error %||% ""), 1, 160),
              proc.time()[["elapsed"]] - t0))
  invisible(r)
}

# 1. the ordinal mixture's centering shift
with_x <- function(s) {
  s$bf <- "bf(y ~ 1, mu1 ~ x, mu2 ~ x)"
  s
}
run("ordmix x: control", "p131", spec_edit = with_x)
run("ordmix x: shift sign flipped", "p131", spec_edit = with_x,
    mutate = function(e) {
      v <- e$eq_rules$ord_mix_intercept$value
      e$eq_rules$ord_mix_intercept$value <- function(nm, fit, sdat, ...) {
        out <- v(nm, fit, sdat, ...)
        k <- sub("^Intercept_mu", "", nm)
        x <- sdat[[paste0("X_mu", k)]]
        fe <- e$brms_fe_of(fit, paste0("mu", k))
        cn <- colnames(x)
        shift <- sum(colMeans(x)[cn] * unname(fe[e$brms_coef_to_frm(cn)]))
        cat("   shift for", nm, "=", shift, "\n")
        out$value <- out$value + 2 * shift
        out
      }
    })

# 2. two cs() columns per component
with_cs2 <- function(s) {
  s$bf <- "bf(y ~ 1, mu1 ~ cs(x) + cs(z), mu2 ~ cs(x) + cs(z))"
  s
}
run("ordmix cs2: control", "p131", spec_edit = with_cs2)
run("ordmix cs2: bcs read by column", "p131", spec_edit = with_cs2,
    mutate = function(e) {
      v <- e$eq_rules$bcs_prefixed$value
      e$eq_rules$bcs_prefixed$value <- function(...) {
        out <- v(...)
        m <- out$value
        out$value <- matrix(as.numeric(t(m)), nrow = nrow(m), ncol = ncol(m))
        out
      }
    })

# 3. REML with a 2-D smooth: a bs map that need not be the identity
with_2d <- function(s) {
  s$bf <- "bf(y ~ s(x, z))"
  s
}
run("REML s(x, z): control", "p044", spec_edit = with_2d,
    mutate = function(e) {
      v <- e$eq_bs_logdet
      e$eq_bs_logdet <- function(fit, sdat) {
        ld <- v(fit, sdat)
        cat("   bs_logdet =", format(ld, digits = 10), " Ks =",
            NCOL(sdat[["Xs"]]), "\n")
        ld
      }
    })
run("REML s(x, z): log-det dropped", "p044", spec_edit = with_2d,
    mutate = function(e) {
      e$eq_bs_logdet <- function(fit, sdat) 0
    })
