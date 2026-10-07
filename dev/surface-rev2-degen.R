# Reviewer re-check, B2: pp_mixture() on degenerate two-gaussian mixture
# fits of unimodal data (n = 40, data seeds 1..30), where a component is
# near empty and the standard errors are large. Traces the largest
# logit SD logitnormal_sd() receives, and catches any pp_mixture() error.
#   Rscript dev/surface-rev2-degen.R > dev/surface-rev-out/rev2-degen.txt
source("dev/surface-rev-env.R"); rev_env("lane")
suppressPackageStartupMessages(library(frmtmb))
mx <- 0
trace("logitnormal_sd", where = asNamespace("frmtmb"), print = FALSE,
      tracer = quote(assign("mx", max(c(get("mx", globalenv()),
                                        abs(lse[is.finite(lse)]))),
                            globalenv())))
res <- character()
for (s in 1:30) {
  set.seed(s)
  d <- data.frame(y = rnorm(40))
  f <- tryCatch(suppressWarnings(suppressMessages(
    frm(bf(y ~ 1), family = mixture(gaussian(), gaussian()), data = d))),
    error = function(e) NULL)
  if (is.null(f)) { res <- c(res, "fitfail"); next }
  mx <- 0
  pm <- tryCatch(suppressWarnings(pp_mixture(f)), error = function(e) e)
  if (inherits(pm, "error")) {
    cat("seed", s, "pp_mixture ERROR:", conditionMessage(pm), "| max lse",
        mx, "\n")
    res <- c(res, "pperr"); next
  }
  cat(sprintf("seed %2d: max logit SD %-10.4g | NA Est.Error %d | min p %.3g\n",
              s, mx, sum(is.na(pm[, "Est.Error", ])),
              min(pm[, "Estimate", ])))
  res <- c(res, "ok")
}
print(table(res))
