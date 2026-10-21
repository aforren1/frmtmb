# Which method a loo_*() call on draws reaches, by load order, on a
# chosen build: brms attached after frmtmb.sample, brms not loaded,
# rstantools attached first. The before arm is rellib-r10.
#   Rscript dev/loodraws-probe-dispatch.R <lane|base> <order>
# order: brms | none | rstantools
# Logs: dev/loodraws-log/dispatch-<arm>-<order>.txt
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(if (a[1] == "lane") "C:/Users/adf44/source/r/wt-loodraws-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (a[2] == "rstantools") suppressMessages(library(rstantools))
suppressMessages(library(frmtmb.sample))
if (a[2] == "brms") suppressMessages(library(brms))
cat("lib:", dirname(system.file(package = "frmtmb.sample")),
    format(packageVersion("frmtmb.sample")), "order:", a[2],
    "brms loaded:", isNamespaceLoaded("brms"), "\n")
set.seed(1)
dd <- data.frame(x = rnorm(40))
dd$y <- rnorm(40, 1 + 0.5 * dd$x)
fit <- frmtmb::frm(frmtmb::bf(y ~ x), family = gaussian(), data = dd)
tpl <- fit$frame$par_template
est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
M <- matrix(rep(est, each = 100) + rnorm(100 * length(est), 0, 0.05), 100,
            dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
ds <- structure(list(stanfit = NULL, draws = M, fit = fit),
                class = "frmtmb_draws")
for (g in c("loo_epred", "loo_linpred", "loo_predict",
            "loo_predictive_interval", "loo_R2")) {
  f <- get0(g, mode = "function")
  where <- if (is.null(f)) "none" else
    environmentName(topenv(environment(f)))
  r <- tryCatch({
    v <- suppressWarnings(suppressMessages(f(ds)))
    paste("ok, dim", paste(dim(v), collapse = " x "))
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-24s generic from %-14s %s\n", g, where, substr(r, 1, 90)))
}
g <- get0("loo_model_weights", mode = "function")
r <- tryCatch({
  v <- suppressWarnings(suppressMessages(g(ds, ds)))
  paste("ok,", paste(names(v), collapse = ","))
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat(sprintf("%-24s generic from %-14s %s\n", "loo_model_weights",
            if (is.null(g)) "none" else
              environmentName(topenv(environment(g))), substr(r, 1, 90)))
