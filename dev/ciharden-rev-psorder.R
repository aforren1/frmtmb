# Reviewer: posterior_samples() deprecation warnings by load order.
#   Rscript dev/ciharden-rev-psorder.R <lib|base> <scenario>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
sc <- a[2]
pre <- switch(sc,
  none = character(),
  brms_ns_after = character(),
  brms_attached_first = "brms",
  brms_attached_after = character(),
  gratia_and_brms = character(),
  stop("scenario"))
for (p in pre) suppressPackageStartupMessages(library(p, character.only = TRUE))
suppressPackageStartupMessages(library(frmtmb.sample))
if (sc == "brms_ns_after") loadNamespace("brms")
if (sc == "brms_attached_after") {
  suppressPackageStartupMessages(library(brms))
}
if (sc == "gratia_and_brms") {
  loadNamespace("gratia")
  loadNamespace("brms")
}
set.seed(20261006)
dd <- data.frame(x = rnorm(30))
dd$y <- rnorm(30, 1 + 0.5 * dd$x)
fit <- frm(frmtmb::bf(y ~ x), family = gaussian(), data = dd)
tpl <- fit$frame$par_template
est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
M <- matrix(rep(est, each = 6) + stats::rnorm(6 * length(est), 0, 0.05),
            6, dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
ds <- structure(list(stanfit = NULL,
                     draws = cbind(frmtmb.sample:::draws_to_natural(M, fit),
                                   lp__ = 0),
                     fit = fit), class = "frmtmb_draws")
count <- function(expr) {
  n <- 0L
  msgs <- character()
  withCallingHandlers(expr, warning = function(w) {
    n <<- n + 1L
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  paste0(n, if (n) paste0(" [", paste(substr(msgs, 1, 40), collapse = " | "),
                          "]"))
}
owner <- function(f) {
  e <- environment(f)
  if (is.null(e)) "?" else environmentName(e)
}
cat("scenario", sc, "| brms loaded:", isNamespaceLoaded("brms"),
    "| gratia loaded:", isNamespaceLoaded("gratia"), "\n")
cat("  search-path posterior_samples from:", owner(posterior_samples), "\n")
cat("  posterior_samples(ds):", count(posterior_samples(ds)), "\n")
cat("  frmtmb.sample::posterior_samples(ds):",
    count(frmtmb.sample::posterior_samples(ds)), "\n")
cat("  lapply(list(ds), posterior_samples):",
    count(lapply(list(ds), posterior_samples)), "\n")
cat("  do.call:", count(do.call("posterior_samples", list(ds))), "\n")
cat("  method direct:",
    count(frmtmb.sample:::posterior_samples.frmtmb_draws(ds)), "\n")
if (isNamespaceLoaded("brms")) {
  cat("  brms::posterior_samples(ds):", count(brms::posterior_samples(ds)),
      "\n")
}
if (isNamespaceLoaded("gratia")) {
  cat("  gratia::posterior_samples(ds):",
      count(tryCatch(gratia::posterior_samples(ds),
                     error = function(e) cat("error:", conditionMessage(e)))),
      "\n")
}
