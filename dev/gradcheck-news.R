# Prepend the development-version heading and this lane's bullet to
# NEWS.md, in one place, so the file is edited by a script rather than by
# a shell quoting dance.

p <- "NEWS.md"
old <- readLines(p, encoding = "UTF-8", warn = FALSE)
if (any(grepl("^# frmtmb \\(development version\\)", old))) {
  stop("the development heading is already there")
}
new <- c(
"# frmtmb (development version)",
"",
"- The convergence warning \"Large maximum absolute gradient at the",
"  optimum\" no longer fires on correct fits. It used to read one",
"  absolute number, `frmtmb_control(grad_tol = 1e-3)`, against the",
"  largest gradient component, and that number does not mean the same",
"  thing on every design: the optimizer stops on a RELATIVE change in an",
"  objective that grows with the sample size, so the gradient at a good",
"  optimum grows with it too, and a parameter held on a bound by",
"  `set_prior(ub = )` has a nonzero gradient by construction. Over 720",
"  correct fits on nine designs at four sample sizes the old check",
"  warned on 289 of them, rising from 11 percent of the fits at 200 rows",
"  to 75 percent at 20,000. It now warns on none of them. `grad_tol` is",
"  read twice: first as the gradient trip-wire it always was, then on",
"  the objective's own scale. A fit that trips the trip-wire warns only",
"  if one exact Newton step over the parameters no bound holds would",
"  still gain more than `grad_tol` in log-likelihood, and the warning",
"  reports that number, so it says how far short the fit is rather than",
"  only that a derivative is large. Every fit's estimates, objective and",
"  optimizer status are bitwise unchanged; only the warning and",
"  `diagnose()` differ. `diagnose()` gains `grad_proj`,",
"  `grad_bound_held` and `grad_headroom`, names the bounds that hold a",
"  gradient, and judges the gradient by `grad_tol` instead of a",
"  hardcoded `1e-3` that ignored the setting. The new criterion also",
"  reports what the gradient could not: a near-collinear design that",
"  nlminb calls converged at a gradient of `5.1e-3` is short by 0.528",
"  log-likelihood units, which the Newton step confirms to 0.2 percent",
"  and which a re-optimization from the same point does not recover.",
"")
writeLines(c(new, old), p, useBytes = TRUE)
cat("wrote", length(new), "lines above", old[[1L]], "\n")
