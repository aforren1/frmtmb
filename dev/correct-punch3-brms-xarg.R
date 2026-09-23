# Punch round 1, minor 1 of the recheck: what does brms do with `x` on a
# bayesplot type that has no `x` argument? The `group` half was measured
# and copied; this half was not measured at all, and a rule read off
# brms's source is not a measurement.
#
# brms's pp_check.brmsfit() builds its bayesplot call with
#   if (!is.null(group)) ppc_args$group <- data[[group]]
#   if (!is.null(x))     ppc_args$x     <- as.numeric(data[[x]])
# so a name the data does not carry gives NULL for `group`, which drops
# out of the list, and `as.numeric(NULL)` for `x`, which is numeric(0)
# and STAYS. That asymmetry is what this measures rather than asserts.
#
#   Rscript dev/correct-punch3-brms-xarg.R
source("C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-prelude.R")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), " bayesplot",
    format(packageVersion("bayesplot")), "\n")
FITS <- "C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-brmsfits"
f <- readRDS(file.path(FITS, "gauss-g.rds"))
cat("fit:", deparse1(formula(f)$formula), " data cols:",
    paste(names(f$data), collapse = " "), "\n\n")

show <- function(label, expr) {
  r <- tryCatch({
    v <- withCallingHandlers(force(expr), warning = function(w) {
      cat(sprintf("%-46s WARNING: %s\n", label,
                  substr(conditionMessage(w), 1, 80)))
      invokeRestart("muffleWarning")
    })
    if (inherits(v, "ggplot")) invisible(ggplot2::ggplot_build(v))
    "OK"
  }, error = function(e) {
    paste0("ERROR: ", substr(gsub("\n", " ", conditionMessage(e)), 1, 110))
  })
  cat(sprintf("%-46s %s\n", label, r))
}

# dens_overlay has neither a `group` nor an `x` formal, so both names
# reach bayesplot's dots or do not reach it at all
cat("== dens_overlay: no group formal, no x formal ==\n")
fa <- names(formals(bayesplot::ppc_dens_overlay))
cat("ppc_dens_overlay formals:", paste(fa, collapse = " "), "\n")
show("group = g      (name IS in the data)", suppressMessages(
  pp_check(f, type = "dens_overlay", ndraws = 5, group = "g")))
show("group = nosuch (name is NOT in the data)", suppressMessages(
  pp_check(f, type = "dens_overlay", ndraws = 5, group = "nosuch")))
show("x = x          (name IS in the data)", suppressMessages(
  pp_check(f, type = "dens_overlay", ndraws = 5, x = "x")))
show("x = nosuch     (name is NOT in the data)", suppressMessages(
  pp_check(f, type = "dens_overlay", ndraws = 5, x = "nosuch")))

# and the same two names on a type that DOES take them, as a control:
# there brms looks the name up and refuses one the data does not carry
cat("\n== controls: a type that takes the argument ==\n")
show("stat_grouped, group = nosuch", suppressMessages(
  pp_check(f, type = "stat_grouped", ndraws = 5, group = "nosuch")))
show("intervals, x = nosuch", suppressMessages(
  pp_check(f, type = "intervals", ndraws = 5, x = "nosuch")))
