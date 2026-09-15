# REVIEW of lane eamhier, ranked item 4: the two-route consistency
# expectation added to tests/testthat/test-scale.R.
#
#   expect_lt(abs(ndt_hat - mean(ndt_sub)),
#             3 * (stats::sd(ndt_sub) + ndt_se))
#
# The lane SAW it fail on the ungrouped row, which is the row the bug
# was on. This script constructs the case where the guarded thing is
# ABSENT on the GROUPED rows as well, because a guard that only has
# power on the row that was already broken is a guard with one seed in
# it. Each mutation below is a way the same class of error can be
# written; the guard has to fire on each.
#
# Run: Rscript --vanilla dev/eamhierrev-guard.R

.libPaths(c("C:/Users/adf44/source/r/eamhierrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(testthat)
  library(frmtmb)
  library(frmtmb.eam)
})
Sys.setenv(NOT_CRAN = "true")
here <- getwd()
setwd("extensions/frmtmb.eam/tests/testthat")
source("helper-scale.R")
# the tier's own data and formula constructors, taken from the file
# under review rather than retyped
src <- readLines("test-scale.R")
i0 <- grep("^eam_truth <- list", src)[1L]
i1 <- grep("^eam_scale_run <- function", src)[1L] - 1L
eval(parse(text = paste(src[i0:i1], collapse = "\n")))
setwd(here)

run_one <- function(group, small = TRUE, sv = 0) {
  Sys.setenv(FRMTMB_SCALE_SMALL = if (small) "true" else "false")
  d <- eam_scale_data(sv = sv)
  form <- eam_scale_form(group)
  fam <- if (sv > 0) wiener(variability = "sv") else wiener()
  t0 <- Sys.time()
  fit <- frm(form, family = fam, data = d, se = TRUE)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  nd <- suppressWarnings(
    stats::predict(fit, newdata = d[1L, , drop = FALSE], dpar = "ndt",
                   type = "response", re.form = NA, se.fit = TRUE))
  bnd <- frmtmb::single_response(fit)[["family"]][["ndt_bound"]]
  one <- d[match(levels(d$s), as.character(d$s)), , drop = FALSE]
  list(d = d, fit = fit, el = el,
       frac = as.numeric(nd$fit[1L]),
       frac_se = as.numeric(nd$se.fit[1L]),
       floors = bnd[["floors"]], ub = bnd[["ub"]],
       ndt_sub = suppressWarnings(as.numeric(ndt_time(fit,
                                                      newdata = one))))
}

# the guard, exactly as the file under review writes it
guard <- function(ndt_hat, ndt_se, ndt_sub) {
  lhs <- abs(ndt_hat - mean(ndt_sub))
  rhs <- 3 * (stats::sd(ndt_sub) + ndt_se)
  c(lhs = lhs, rhs = rhs, fires = as.numeric(lhs >= rhs),
    margin = rhs / lhs)
}

report <- function(tag, r) {
  cat("\n== ", tag, " ==  ", sprintf("%.1f s", r$el), "\n", sep = "")
  cat(sprintf("  subjects %d | floors %s | ub %s\n",
              length(r$ndt_sub),
              if (is.null(r$floors)) "NULL" else
                sprintf("%d, mean %.6f", length(r$floors),
                        mean(r$floors)),
              if (is.null(r$ub)) "NULL" else sprintf("%.6f", r$ub)))
  cat(sprintf(paste0("  predict() fraction/time %.6f (se %.6f) ",
                     "| mean ndt_time() %.6f | sd %.6f\n"),
              r$frac, r$frac_se, mean(r$ndt_sub),
              stats::sd(r$ndt_sub)))
  to_time <- if (is.null(r$floors)) 1 else mean(r$floors)
  muts <- list(
    "as shipped (fixed)" = to_time,
    "PRE-FIX idiom (ub when ungrouped)" =
      if (is.null(r$floors)) r$ub else mean(r$floors),
    "forgot the conversion (to_time = 1)" = 1,
    "converted twice (to_time^2)" = to_time^2)
  for (nm in names(muts)) {
    tt <- muts[[nm]]
    g <- guard(r$frac * tt, r$frac_se * tt, r$ndt_sub)
    cat(sprintf(paste0("  %-36s to_time %8.6f  |diff| %9.6f  ",
                       "tol %9.6f  %s (tol/|diff| %6.2f)\n"),
                nm, tt, g[["lhs"]], g[["rhs"]],
                if (g[["fires"]] > 0) "FIRES " else "passes",
                g[["margin"]]))
  }
}

report("SMALL, grouped (rows eam / eam-sv)", run_one(TRUE))
report("SMALL, ungrouped (row eam-unbounded)", run_one(FALSE))

# The full-size grouped row, from the numbers this lane recorded for it
# rather than from a fresh 110 s fit: seed 20260908, tier row `eam`.
cat("\n== FULL SIZE, grouped, from dev/eamhier-findings.md's own record ==\n")
frac <- 0.2468670903 / 0.290527
sub_mean <- 0.2461460561
sub_sd <- 0.025425
se <- 0.0019639
for (nm in c("as shipped", "to_time = 1", "converted twice")) {
  tt <- switch(nm, "as shipped" = 0.290527, "to_time = 1" = 1,
               "converted twice" = 0.290527^2)
  lhs <- abs(frac * tt - sub_mean)
  rhs <- 3 * (sub_sd + frac * tt * (se / (frac * 0.290527)))
  cat(sprintf("  %-18s |diff| %9.6f  tol %9.6f  %s\n", nm, lhs, rhs,
              if (lhs >= rhs) "FIRES" else "passes"))
}
cat("\ndone\n")
