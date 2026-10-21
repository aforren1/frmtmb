# Reviewer of lane loodraws, fix round 1: the draws-only names left in
# frmtmb.sample (loo_moment_match, loo_subsample) and the moved kfold and
# reloo, on an ML fit, by attach order; and row 469's `fit6$data`.
#   Rscript dev/loodraws-rev-gap.R <rev|base> <order>
# order: frmtmb | loo-frmtmb | brms-frmtmb | sample
# Log: dev/loodraws-rev-out/gap-<arm>-<order>.txt
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(if (a[1] == "rev") "C:/Users/adf44/source/r/wt-loodraws-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
att <- function(p) suppressMessages(library(p, character.only = TRUE))
switch(a[2],
  frmtmb = att("frmtmb"),
  `loo-frmtmb` = { att("loo"); att("frmtmb") },
  `brms-frmtmb` = { att("brms"); att("frmtmb") },
  sample = att("frmtmb.sample"),
  `rst-frmtmb` = { att("rstantools"); att("frmtmb") },
  `sample-loo` = { att("frmtmb.sample"); att("loo") },
  `loo-sample` = { att("loo"); att("frmtmb.sample") },
  `sample-brms` = { att("frmtmb.sample"); att("brms") },
  `brms-sample` = { att("brms"); att("frmtmb.sample") },
  loadedonly = { att("frmtmb"); invisible(loadNamespace("frmtmb.sample")) },
  stop("unknown order"))
cat("arm", a[1], "order", a[2], "| search:",
    paste(grep("^package:(frmtmb|brms|loo|rstantools)", search(),
               value = TRUE), collapse = " "), "\n")
set.seed(1)
dd <- data.frame(x = rnorm(40))
dd$y <- rnorm(40, 1 + 0.5 * dd$x)
fit <- suppressWarnings(frmtmb::frm(frmtmb::bf(y ~ x), family = gaussian(),
                                    data = dd))
for (g in c("kfold", "reloo", "loo_moment_match", "loo_subsample")) {
  f <- get0(g, mode = "function")
  where <- if (is.null(f)) "none" else environmentName(topenv(environment(f)))
  r <- if (is.null(f)) "could not find function" else tryCatch({
    suppressWarnings(suppressMessages(f(fit))); "ok"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("fit %-18s from %-12s %s\n", g, where,
              substr(gsub("\n", " ", r), 1, 130)))
}
if (isNamespaceLoaded("frmtmb.sample")) {
  tpl <- fit$frame$par_template
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  M <- matrix(rep(est, each = 100) + rnorm(100 * length(est), 0, 0.05),
              100, dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  dz <- structure(list(stanfit = NULL, draws = M, fit = fit),
                  class = "frmtmb_draws")
  for (g in c("loo_moment_match", "loo_subsample", "kfold", "reloo")) {
    fs <- c(list(bare = get0(g, mode = "function")),
            if (g %in% getNamespaceExports("loo"))
              list(`loo::` = getExportedValue("loo", g)),
            if (g %in% getNamespaceExports("brms"))
              list(`brms::` = getExportedValue("brms", g)))
    for (nm in names(fs)) {
      f <- fs[[nm]]
      r <- if (is.null(f)) "could not find function" else tryCatch({
        suppressWarnings(suppressMessages(f(dz))); "ok"
      }, error = function(e) paste("ERROR:", conditionMessage(e)))
      cat(sprintf("draws %-6s %-18s %s\n", nm, g,
                  substr(gsub("\n", " ", r), 1, 120)))
    }
  }
}
if (a[2] == "sample") {
  tpl <- fit$frame$par_template
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  M <- matrix(rep(est, each = 100) + rnorm(100 * length(est), 0, 0.05),
              100, dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  ds <- structure(list(stanfit = NULL, draws = M, fit = fit),
                  class = "frmtmb_draws")
  cat("row 469: names(ds):", names(ds), "| is.null(ds$data):",
      is.null(ds$data), "\n")
  e1 <- suppressWarnings(suppressMessages(loo(ds, cores = 1)))
  e2 <- suppressWarnings(suppressMessages(loo(ds, cores = 1,
                                              newdata = ds$data)))
  cat("row 469: loo(newdata = ds$data) identical to loo():",
      identical(e1$estimates, e2$estimates), "\n")
  r <- tryCatch({
    suppressWarnings(suppressMessages(loo(ds, newdata = model.frame(ds))))
    "ok"
  }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 90)))
  cat("row 469 with the data given: loo(newdata = model.frame(ds)):", r, "\n")
}
