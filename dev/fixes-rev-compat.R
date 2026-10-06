# Reviewer of lane fixes, claim 4d: the compat cells the lane's probe
# did not run (disc intercept warning, delta prior, sum_to_zero
# tau_raw prior, equidistant coef refusal).
#   Rscript dev/fixes-rev-compat.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
set.seed(62)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n))
u <- stats::rlogis(n) / exp(0.3 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
tr <- function(lab, expr) {
  w <- character()
  r <- tryCatch(withCallingHandlers({expr; "OK"}, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 150)))
  cat(sprintf("%-46s %s | warnings: %s\n", lab, r,
              if (length(w)) substr(paste(w, collapse = " / "), 1, 150) else "none"))
}
tr("disc ~ 1 + z (intercept in disc)",
   frm(bf(y ~ x, disc ~ 1 + z), family = cumulative(), data = d))
tr("equidistant, class delta prior",
   frm(y ~ x, family = cumulative(threshold = "equidistant"), data = d,
       prior = set_prior("normal(1, 0.5)", class = "delta", lb = 0)))
tr("equidistant, class Intercept prior",
   frm(y ~ x, family = cumulative(threshold = "equidistant"), data = d,
       prior = set_prior("normal(0, 1)", class = "Intercept")))
tr("sum_to_zero, tau_raw prior",
   frm(y ~ x, family = sratio(threshold = "sum_to_zero"), data = d,
       prior = list(tau_raw = prior_normal(0, 2))))
tr("sum_to_zero, class Intercept prior",
   frm(y ~ x, family = sratio(threshold = "sum_to_zero"), data = d,
       prior = set_prior("normal(0, 1)", class = "Intercept")))
tr("disc prior, dpar = disc class b",
   frm(bf(y ~ x, disc ~ 0 + z), family = cumulative(), data = d,
       prior = set_prior("normal(0, 1)", class = "b", dpar = "disc")))
for (st in c("disc", "equidistant", "sum_to_zero")) {
  r <- frm_compat(st)
  cat(st, ": ", paste(r$with, r$status, sep = "=", collapse = " "), "\n",
      sep = "")
}
