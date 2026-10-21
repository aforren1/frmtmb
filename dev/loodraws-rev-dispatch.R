# Reviewer of lane loodraws: what a user sees from the loo_*() names, on a
# draws object and on an ML fit, by what is attached and in what order.
#   Rscript dev/loodraws-rev-dispatch.R <rev|base> <order>
# order (left to right is attach order):
#   frmtmb        frmtmb alone; frmtmb.sample never loaded
#   rst-frmtmb    rstantools, then frmtmb alone
#   brms-frmtmb   brms, then frmtmb alone
#   sample-loo    frmtmb.sample, then loo
#   loo-sample    loo, then frmtmb.sample
#   sample-rst    frmtmb.sample, then rstantools
#   loadedonly    frmtmb.sample loaded, not attached; frmtmb attached
# Log: dev/loodraws-rev-out/dispatch-<arm>-<order>.txt
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(if (a[1] == "rev") "C:/Users/adf44/source/r/wt-loodraws-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ord <- a[2]
att <- function(p) suppressMessages(library(p, character.only = TRUE))
switch(ord,
  frmtmb = att("frmtmb"),
  `rst-frmtmb` = { att("rstantools"); att("frmtmb") },
  `brms-frmtmb` = { att("brms"); att("frmtmb") },
  `sample-loo` = { att("frmtmb.sample"); att("loo") },
  `loo-sample` = { att("loo"); att("frmtmb.sample") },
  `sample-rst` = { att("frmtmb.sample"); att("rstantools") },
  loadedonly = { att("frmtmb"); invisible(loadNamespace("frmtmb.sample")) },
  `sample-brms` = { att("frmtmb.sample"); att("brms") },
  `brms-sample` = { att("brms"); att("frmtmb.sample") },
  `loo-frmtmb` = { att("loo"); att("frmtmb") },
  stop("unknown order ", ord))
cat("arm", a[1], "order", ord, "| sample loaded",
    isNamespaceLoaded("frmtmb.sample"), "| brms loaded",
    isNamespaceLoaded("brms"), "| search:",
    paste(grep("^package:(frmtmb|brms|loo|rstantools)", search(),
               value = TRUE), collapse = " "), "\n")
set.seed(1)
dd <- data.frame(x = rnorm(40))
dd$y <- rnorm(40, 1 + 0.5 * dd$x)
fit <- frmtmb::frm(frmtmb::bf(y ~ x), family = gaussian(), data = dd)
ds <- NULL
if (isNamespaceLoaded("frmtmb.sample")) {
  tpl <- fit$frame$par_template
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  M <- matrix(rep(est, each = 100) + rnorm(100 * length(est), 0, 0.05),
              100, dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  ds <- structure(list(stanfit = NULL, draws = M, fit = fit),
                  class = "frmtmb_draws")
}
show <- function(g, obj, label) {
  f <- get0(g, mode = "function")
  where <- if (is.null(f)) "none" else environmentName(topenv(environment(f)))
  r <- if (is.null(f)) "could not find function" else tryCatch({
    v <- suppressWarnings(suppressMessages(
      if (g == "loo_model_weights") f(obj, obj) else f(obj)))
    paste("ok", paste(dim(v), collapse = "x"))
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-6s %-24s from %-14s %s\n", label, g, where,
              substr(gsub("\n", " ", r), 1, 150)))
}
for (g in c("loo_epred", "loo_linpred", "loo_predict",
            "loo_predictive_interval", "loo_R2", "loo_model_weights",
            "loo", "bayes_R2")) {
  show(g, fit, "fit")
  if (!is.null(ds)) show(g, ds, "draws")
}
# the owners' generics called by their full names
for (pk in c("rstantools", "brms", "loo")) {
  for (g in c("loo_epred", "loo_R2", "loo_model_weights")) {
    if (!g %in% getNamespaceExports(pk)) next
    f <- getExportedValue(pk, g)
    for (o in list(list(fit, "fit"), list(ds, "draws"))) {
      if (is.null(o[[1]])) next
      r <- tryCatch({
        v <- suppressWarnings(suppressMessages(
          if (g == "loo_model_weights") f(o[[1]], o[[1]]) else f(o[[1]])))
        paste("ok", paste(dim(v), collapse = "x"))
      }, error = function(e) paste("ERROR:", conditionMessage(e)))
      cat(sprintf("%-6s %s::%-20s %s\n", o[[2]], pk, g,
                  substr(gsub("\n", " ", r), 1, 150)))
    }
  }
}
