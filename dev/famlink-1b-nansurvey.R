# Round 1b, item 2: which (family, link) pairs brms accepts raise
# "NA/NaN function evaluation" on a correct fit, and does the fix change
# any optimum. Same pairs, data, seed and formula as
# dev/famlink-falsealarm.R (seed 20260916, n = 120), whose data and
# constructor code this script sources rather than copies.
#
#   Rscript dev/famlink-1b-nansurvey.R before   (lane build before the fix)
#   Rscript dev/famlink-1b-nansurvey.R after    (lane build after the fix)
#
# writes dev/famlink-1b-nansurvey-<tag>.rds: per pair the NaN warning
# count, the other warnings, logLik and the estimate vector.
tag <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(tag %in% c("before", "after"))
.libPaths(c("C:/Users/adf44/source/r/famlink-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))

# the data builders of the false-alarm script, taken from its source
src <- parse("dev/famlink-falsealarm.R", keep.source = FALSE)
keep <- vapply(src, function(e) {
  is.call(e) && identical(e[[1]], as.name("<-")) &&
    as.character(e[[2]]) %in% c("n", "make_data", "form_for", "ctor_for")
}, NA)
for (e in src[keep]) eval(e, globalenv())
tab_env <- new.env()
sys.source("R/links-brms.R", envir = tab_env)
sets <- tab_env$brms_mu_links

out <- list()
for (fam in setdiff(names(sets), "multinomial")) {
  d <- make_data(fam)
  for (lk in sets[[fam]]) {
    warns <- character(0)
    r <- withCallingHandlers(tryCatch({
      f <- frmtmb:::family_registry[[ctor_for(fam)]](link = lk)
      fit <- suppressMessages(frm(form_for(fam), d, family = f))
      list(ll = as.numeric(logLik(fit)), par = fit$opt$par,
           conv = fit$opt$convergence)
    }, error = function(e) list(error = conditionMessage(e))),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    nan <- sum(grepl("NA/NaN function evaluation", warns, fixed = TRUE))
    out[[paste(fam, lk)]] <- c(r, list(family = fam, link = lk, nan = nan,
                                       other = setdiff(warns,
      "NA/NaN function evaluation")))
    cat(sprintf("%-28s %-14s NaN warnings %2d %s\n", fam, lk, nan,
                if (!is.null(r$error)) "ERROR" else ""))
  }
}
saveRDS(out, paste0("dev/famlink-1b-nansurvey-", tag, ".rds"))
cat("pairs:", length(out), " pairs with a NaN warning:",
    sum(vapply(out, function(o) o$nan > 0, NA)), "\n")

if (tag == "after" && file.exists("dev/famlink-1b-nansurvey-before.rds")) {
  bef <- readRDS("dev/famlink-1b-nansurvey-before.rds")
  stopifnot(identical(names(bef), names(out)))
  same <- vapply(names(out), function(k) {
    identical(bef[[k]]$ll, out[[k]]$ll) && identical(bef[[k]]$par, out[[k]]$par)
  }, NA)
  cat("---- GENERATED: dev/famlink-1b-nansurvey.R ----\n")
  cat("pairs:", length(out), "\n")
  cat("with a NaN warning before:", sum(vapply(bef, function(o) o$nan > 0, NA)),
      "; after:", sum(vapply(out, function(o) o$nan > 0, NA)), "\n")
  cat("logLik and every outer estimate identical() before and after:",
      sum(same), "of", length(same), "\n")
  cat("---- END GENERATED ----\n")
}
