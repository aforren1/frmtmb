# Reviewer: compare warmup-old.rds and warmup-new.rds.
L <- "C:/Users/adf44/source/r/frmtmb-wt-priordraws/dev/review-priordraws-log"
o <- readRDS(file.path(L, "warmup-old.rds")); n <- readRDS(file.path(L, "warmup-new.rds"))
for (nm in names(o)) {
  cat(sprintf("%-10s draws identical %s | inc_warmup array identical %s | rhat(inc_warmup) identical %s | as.mcmc dims identical %s",
              nm, identical(o[[nm]]$draws, n[[nm]]$draws), identical(o[[nm]]$a, n[[nm]]$a),
              identical(o[[nm]]$rh, n[[nm]]$rh), identical(o[[nm]]$mc, n[[nm]]$mc)))
  if (!is.null(o[[nm]]$size)) cat(sprintf(" | size %.0f -> %.0f KB", o[[nm]]$size / 1024, n[[nm]]$size / 1024))
  if (!identical(o[[nm]]$a, n[[nm]]$a) && is.array(o[[nm]]$a) && is.array(n[[nm]]$a)) {
    w <- seq_len(dim(o[[nm]]$a)[1] - nrow(o[[nm]]$draws) / dim(o[[nm]]$a)[2])
    cat(sprintf(" | max abs diff warmup rows %.3g", max(abs(o[[nm]]$a[w, , ] - n[[nm]]$a[w, , ]))))
  }
  cat("\n")
}
