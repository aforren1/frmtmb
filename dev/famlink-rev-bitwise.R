## Reviewer check for lane wt-famlink, priority 4: the worker's 207
## "logLik equal to 10 printed decimals" fits, compared BITWISE on
## logLik, every estimate, every SE and fitted values.
##
## The data, formula and constructor for each pair are taken from
## dev/famlink-falsealarm.R by parsing out its make_data(), form_for(),
## ctor_for() and `n`, so this is the worker's construction (seed
## 20260916, n = 120) with a stricter comparison.
## Usage: Rscript dev/famlink-rev-bitwise.R <lane|base>
##        Rscript dev/famlink-rev-bitwise.R compare
arm <- commandArgs(trailingOnly = TRUE)[1]
if (identical(arm, "compare")) {
  b <- readRDS("dev/famlink-rev-bitwise-base.rds")
  l <- readRDS("dev/famlink-rev-bitwise-lane.rds")
  both <- intersect(names(b)[vapply(b, function(r) r$status == "fit", TRUE)],
                    names(l)[vapply(l, function(r) r$status == "fit", TRUE)])
  flds <- c("ll", "est", "se", "fitted")
  nd <- 0L
  for (k in both) {
    same <- vapply(flds, function(f) identical(b[[k]][[f]], l[[k]][[f]]), TRUE)
    if (!all(same)) {
      nd <- nd + 1L
      cat(k, "differs in", flds[!same], " max |d est|",
          max(abs(b[[k]]$est - l[[k]]$est)), "\n")
    }
  }
  cat("pairs:", length(b), " fitted in both:", length(both),
      " bitwise identical in logLik, estimates, SEs, fitted:",
      length(both) - nd, "\n")
  cat("fitted base only:", sum(vapply(b, function(r) r$status == "fit", TRUE)) - length(both),
      " fitted lane only:", sum(vapply(l, function(r) r$status == "fit", TRUE)) - length(both), "\n")
  quit(save = "no")
}
ARM <- arm
source("dev/famlink-rev-common.R")
ex <- parse("dev/famlink-falsealarm.R", keep.source = FALSE)
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      as.character(e[[2]]) %in% c("n", "make_data", "form_for", "ctor_for")) {
    eval(e, globalenv())
  }
}
tab_env <- new.env()
sys.source("R/links-brms.R", envir = tab_env)
sets <- tab_env$brms_mu_links
out <- list()
for (fam in setdiff(names(sets), "multinomial")) {
  d <- make_data(fam)
  for (lk in sets[[fam]]) {
    key <- paste(fam, lk)
    out[[key]] <- tryCatch(suppressWarnings(suppressMessages({
      f <- frmtmb:::family_registry[[ctor_for(fam)]](link = lk)
      fit <- frm(form_for(fam), d, family = f)
      list(status = "fit", ll = as.numeric(logLik(fit)),
           est = unlist(fit$estimates),
           se = tryCatch(sqrt(diag(vcov(fit))), error = function(e) NULL),
           fitted = tryCatch(unlist(fitted(fit)), error = function(e) NULL))
    })), error = function(e) list(status = "error"))
  }
}
saveRDS(out, sprintf("dev/famlink-rev-bitwise-%s.rds", arm))
cat("arm", arm, length(out), "pairs\n")
