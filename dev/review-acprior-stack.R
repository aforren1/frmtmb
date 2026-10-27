# Reviewer: keep_flat_defaults() against frm_sample()'s three-source
# stack (call prior, the fit's own prior, the defaults).
# base = fit[["prior"]], as frm_sample() passes it.
# Usage: Rscript dev/review-acprior-stack.R before|after
# Log: dev/review-acprior-log/stack-<arm>.txt
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (identical(arm, "after")) "C:/Users/adf44/source/r/wt-acprior-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(5)
d <- data.frame(t = rep(1:5, 12), g = rep(1:12, each = 5), x = rnorm(60))
d$y <- d$x + rnorm(12)[d$g] * 0.5 + rnorm(60)
fmt <- function(v) {
  if (!length(v)) return("(none)")
  paste(names(v), signif(unlist(v), 6), collapse = "; ")
}
show <- function(lab, fit, pr) {
  cat("\n---", lab, "\n")
  tryCatch({
    m <- character(0)
    ri <- withCallingHandlers(
      frmtmb.sample:::sample_resolve_priors(fit, pr, base = fit[["prior"]])$ri,
      message = function(x) {
        m <<- c(m, conditionMessage(x))
        invokeRestart("muffleMessage")
      })
    for (e in ri$entries) {
      if (!identical(e$comp, "thetaac")) next
      cat("  entry thetaac", paste(e$idx, collapse = ","), e$dist$kind,
          if (identical(e$dist$kind, "trans")) e$dist$inner$kind, "\n")
    }
    cat("  lower:", fmt(ri$lower), " upper:", fmt(ri$upper), "\n")
    ln <- grep("cosy|ar ", unlist(strsplit(m, "\n")), value = TRUE)
    if (length(ln)) cat(paste0("  msg| ", ln), sep = "\n")
  }, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))
}
fc <- frm(y ~ x + cosy(t, g), data = d,
          prior = set_prior("normal(0, 0.1)", class = "cosy"))
show("MAP fit with cosy normal(0, 0.1), no call prior", fc, NULL)
fa <- frm(y ~ x + ar(t, g, cov = TRUE), data = d,
          prior = set_prior("normal(0.5, 0.2)", class = "ar"))
show("MAP fit with ar normal(0.5, 0.2), call prior ar lb = 0.5", fa,
     set_prior("", class = "ar", lb = 0.5))
fb <- frm(y ~ x + ar(t, g, cov = TRUE), data = d,
          prior = set_prior("", class = "ar", lb = 0.2))
show("ML fit with ar lb = 0.2 (bounds only), no call prior", fb, NULL)
show("ML fit with ar lb = 0.2, call prior ar normal(0.5, 0.2)", fb,
     set_prior("normal(0.5, 0.2)", class = "ar"))
