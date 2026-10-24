# Release review of 0.74.0, a cross-lane case no lane ran: a
# mixture_mvn() fit saved by 0.74.0 with frmtmb.latent 0.7.0 (lane
# mvnlatent's move), read where frmtmb.latent is not installed (lane
# splitsurv's require_model_packages()), and a mixture() of a core
# family and a royston_parmar() component (the component-package list).
#   Rscript dev/review-rel074-mvnsaved.R save <dir>
#   Rscript dev/review-rel074-mvnsaved.R read <dir> <lib;..>
a <- commandArgs(trailingOnly = TRUE)
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
if (a[1] == "save") {
  .libPaths(c("C:/Users/adf44/source/r/rellib-r12", user))
  suppressMessages({library(frmtmb); library(frmtmb.latent)})
  dm <- data.frame(row = seq_len(nrow(datasets::faithful)))
  dm$Y <- data.matrix(datasets::faithful)
  fmvn <- frm(bf(Y ~ 1), family = mixture_mvn(K = 2, D = 2), data = dm)
  cat("fmvn logLik", format(as.numeric(logLik(fmvn)), digits = 12),
      "| family package attr:",
      attr(fmvn$spec$responses[[1]]$family, "frmtmb_package", exact = TRUE), "\n")
  saveRDS(list(fit = fmvn, data = dm), file.path(a[2], "fmvn-074.rds"))
  quit(save = "no")
}
.libPaths(c(strsplit(a[3], ";", fixed = TRUE)[[1]], user))
suppressMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), dirname(find.package("frmtmb")),
    "| latent installed:", nzchar(system.file(package = "frmtmb.latent")), "\n")
s <- readRDS(file.path(a[2], "fmvn-074.rds"))
fit <- s$fit
calls <- list(
  print = quote(utils::capture.output(print(fit))),
  summary = quote(utils::capture.output(summary(fit))),
  logLik = quote(logLik(fit)),
  vcov = quote(vcov(fit)),
  fitted = quote(fitted(fit)),
  mixture_probs = quote(mixture_probs(fit)),
  simulate = quote(simulate(fit, nsim = 1, seed = 1)),
  update = quote(update(fit, data = s$data)))
for (k in names(calls)) {
  r <- tryCatch({ eval(calls[[k]]); "ok" }, error = function(e) {
    kk <- if (inherits(e, "frmtmb_extension_missing")) "REFUSED[ext]" else "ERR"
    paste(kk, substr(gsub("\\s+", " ", conditionMessage(e)), 1, 160))
  })
  cat(sprintf("  %-14s %s\n", k, r))
}
