# Lane hollow, Row 2 (brmsfit-methods:469) on the bin-2 draws stand-in
# for brms's example 6, built as helper-brms-suite-draws.R builds it.
#   Rscript dev/hollow-row2-probe.R [base] > dev/hollow-out/row2-<arm>.txt
base <- identical(commandArgs(TRUE)[1], "base")
.libPaths(c(if (!base) "C:/Users/adf44/source/r/wt-hollow-lib",
            "C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-hollow/dev/stan-cache")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), "\n")
cat("lib:", find.package("frmtmb.sample"), "\n")
q <- function(expr) {
  withCallingHandlers(expr, warning = function(w) invokeRestart(
    "muffleWarning"), message = function(m) invokeRestart("muffleMessage"))
}
d <- as.data.frame(brms:::brmsfit_example6$data)
attr(d, "terms") <- NULL
fit <- q(frm(bf(volume ~ Trt + gp(Age, by = Trt, gr = TRUE)) + gaussian() +
               bf(count ~ Trt + Age) + poisson() + set_rescor(FALSE),
             data = d))
fit6 <- q(frm_sample(fit, chains = 1, iter = 75, warmup = 50,
                     seed = 20260917, refresh = 0))
cat("names(fit6):", names(fit6), "\n")
cat("is.null(fit6$data):", is.null(fit6$data), "\n")
loo6_1 <- q(loo(fit6, cores = 1))
loo6_2 <- tryCatch(q(loo(fit6, cores = 1, newdata = fit6$data)),
                   error = function(e) e)
show <- function(label, l) {
  if (inherits(l, "error")) {
    cat(label, ": ERROR ", conditionMessage(l), "\n", sep = "")
  } else {
    cat(label, ": elpd_loo ", format(l$estimates["elpd_loo", 1],
                                     digits = 12),
        " identical to loo(fit6): ",
        identical(l$estimates, loo6_1$estimates), "\n", sep = "")
  }
}
show("loo(fit6)", loo6_1)
show("loo(fit6, newdata = fit6$data)", loo6_2)
show("loo(fit6, newdata = model.frame(fit6))",
     tryCatch(q(loo(fit6, cores = 1, newdata = model.frame(fit6))),
              error = function(e) e))
show("loo(fit6, newdata = d), the data fitted",
     tryCatch(q(loo(fit6, cores = 1, newdata = d)),
              error = function(e) e))
d2 <- d
d2$volume[1] <- d2$volume[1] + 1
show("loo(fit6, newdata = d with one response moved)",
     tryCatch(q(loo(fit6, cores = 1, newdata = d2)),
              error = function(e) e))
cat("ROW2 PROBE END\n")
