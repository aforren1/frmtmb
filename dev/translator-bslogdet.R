# Review 2026-11-09, m3: is there a REML spelling where the map between
# brms's unpenalized smooth columns (Xs) and frmtmb's (<term>.fx<j>) is
# not the identity, so that eq_bs_logdet() is nonzero? For each candidate
# it prints log|det A|, how far the two column sets are apart, and the
# REML marginal identity three ways: as built, with the log-determinant
# dropped, and with its sign flipped. The last candidate is constructed
# so that the log-determinant is not 0 (spec field l1_sdat_fix). Records
# nothing.
#
#   Rscript dev/translator-bslogdet.R
#
# From the worktree root; reads FRMTMB_EQUIV_LIB as dev/equiv-measure.R.
lib <- Sys.getenv("FRMTMB_EQUIV_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE"))) {
  Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache",
                                               winslash = "/"))
}
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

# the survey's data of pair p044 (seed 44), with variants of x
cands <- list(
  "s(x)" = list(bf = "bf(y ~ s(x))"),
  "s(x, m = 3) (two null-space columns)" = list(bf = "bf(y ~ s(x, m = 3))"),
  "t2(x, z)" = list(bf = "bf(y ~ t2(x, z))"),
  "s(x, by = f)" = list(bf = "bf(y ~ f + s(x, by = f))"),
  "s(x) + s(z)" = list(bf = "bf(y ~ s(x) + s(z))"),
  "s(xs), x times 1000 plus 5000" = list(bf = "bf(y ~ s(xs))", big = TRUE),
  "s(x, bs = 'cr')" = list(bf = "bf(y ~ s(x, bs = 'cr'))"),
  # the constructed case: brms's own program read on its own data with
  # the unpenalized columns times 3. That is a change of variables of bs
  # (brms's bs is then frmtmb's / 3), so the joint identity is unchanged
  # and the integral over bs gains log|det A| = -Ks log 3
  "s(x, m = 3), brms's Xs times 3" = list(
    bf = "bf(y ~ s(x, m = 3))",
    sdat_fix = function(sd) {
      sd$Xs <- 3 * sd$Xs
      sd
    })
)
for (nm in names(cands)) {
  cd <- cands[[nm]]
  for (mode in c("asis", "dropped", "flipped")) {
    e <- fresh()
    ld_orig <- e$eq_bs_logdet
    if (mode == "dropped") e$eq_bs_logdet <- function(fit, sdat) 0
    if (mode == "flipped") {
      e$eq_bs_logdet <- function(fit, sdat) -ld_orig(fit, sdat)
    }
    seen <- new.env()
    if (mode == "asis") {
      e$eq_bs_logdet <- function(fit, sdat) {
        v <- ld_orig(fit, sdat)
        seen$ld <- v
        xb <- sdat[["Xs"]]
        if (!is.null(xb)) {
          fe <- e$brms_fe_of(fit, "mu")
          fx <- names(fe)[grepl("[.]fx[0-9]+$", names(fe))]
          xf <- as.matrix(e$brms_lp_of(fit, "mu")$X)[, fx, drop = FALSE]
          seen$gap <- if (identical(dim(as.matrix(xb)), dim(xf))) {
            max(abs(as.matrix(xb) - xf))
          } else NA
          seen$ks <- NCOL(xb)
        }
        v
      }
    }
    d <- e$eq_dat(seed = 44)
    d$y <- sin(2 * d$x) + rnorm(nrow(d), 0, 0.5)
    if (isTRUE(cd$big)) d$xs <- 5000 + 1000 * d$x
    e$eq_specs[["bslogdet x REML"]] <- function() {
      list(bf = cd$bf, data = d, frm_args = list(REML = TRUE),
           l1_marginal = TRUE, l1_sdat_fix = cd$sdat_fix)
    }
    row <- data.frame(id = "m3", feature_a = "bslogdet", feature_b = "REML",
                      status = "untested", stringsAsFactors = FALSE)
    r <- tryCatch(e$eq_measure("bslogdet x REML", row),
                  error = function(err) list(l1_error = conditionMessage(err)))
    cat(sprintf(paste0("%-38s %-7s logdet %s  |Xs - fx| %s  Ks %s  joint %s",
                       "  marginal %s  L1 ok %s %s\n"),
                nm, mode, format(seen$ld %||% NA, digits = 6),
                format(seen$gap %||% NA, digits = 3), seen$ks %||% NA,
                fmt(c(r$l1_const, r$l1_extra_const)), fmt(r$l1_marg),
                format(r$l1_ok %||% NA),
                substr(paste(r$l1_error %||% "", r$frm_error %||% "",
                             r$verdict %||% ""), 1, 120)))
  }
}
