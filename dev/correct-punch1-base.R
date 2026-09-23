# Punch round 1, the BASE arm: what 0.61.0 does with the calls minors 2,
# 6 and 8 changed. The lane's own test file cannot show this, because
# its `pp_check() refuses as brms refuses` block ERRORS on base at the
# `"violin"` assertion, which aborts the block before these three calls
# are reached. A behavioral failure has to be constructed, so it is.
#   Rscript dev/correct-punch1-base.R
source("C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-prelude.R")
.libPaths(.libPaths()[-1L])
suppressMessages(library(frmtmb))
cat("arm base frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n\n")
source("C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-data.R")

show <- function(label, expr) {
  r <- tryCatch({
    v <- withCallingHandlers(force(expr), warning = function(w) {
      cat(sprintf("%-42s WARNING: %s\n", label,
                  substr(conditionMessage(w), 1, 90)))
      invokeRestart("muffleWarning")
    })
    if (inherits(v, "ggplot")) invisible(ggplot2::ggplot_build(v))
    "OK"
  }, error = function(e) {
    paste0("ERROR [", class(e)[1], "]: ",
           substr(gsub("\n", " ", conditionMessage(e)), 1, 120))
  })
  cat(sprintf("%-42s %s\n", label, r))
}

cat("== minor 6: does the refusal quote the caller's type? ==\n")
d <- correct_data_ord()
fo <- frm(bf(ord ~ x) + cumulative(), data = d)
for (ty in c("response", "ordinary", "pearson")) {
  v <- tryCatch(residuals(fo, type = ty),
                error = function(e) conditionMessage(e))
  cat(sprintf("%-12s %s\n", ty,
              if (is.character(v)) substr(v, 1, 70) else
                paste("answers, first 3:",
                      paste(signif(utils::head(v[, 1], 3), 4),
                            collapse = " "))))
}

cat("\n== minor 2: is `resp` a formal at all, and what does it do? ==\n")
cat("formals of pp_check.frmtmb_fit:",
    paste(names(formals(frmtmb:::pp_check.frmtmb_fit)), collapse = " "),
    "\n")
dg <- correct_data_gauss()
fg <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dg)
show("dens_overlay, resp = \"w\"", suppressMessages(
  pp_check(fg, type = "dens_overlay", ndraws = 5, resp = "w")))

cat("\n== minor 8: does bayesplot get the COLUMN, and does it warn? ==\n")
show("dens_overlay, group = g", suppressMessages(
  pp_check(fg, type = "dens_overlay", group = "g", ndraws = 5)))
show("dens_overlay, group = nosuch", suppressMessages(
  pp_check(fg, type = "dens_overlay", group = "nosuch", ndraws = 5)))
