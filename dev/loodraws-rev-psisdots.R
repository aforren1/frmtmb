# Reviewer of lane loodraws: does psis() on draws honor or refuse
# draw_ids and other names in `...`?
#   Rscript dev/loodraws-rev-psisdots.R <rev|base> -> dev/loodraws-rev-out/psisdots-<arm>.txt
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(if (a[1] == "rev") "C:/Users/adf44/source/r/wt-loodraws-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", dirname(system.file(package = "frmtmb.sample")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(1)
dd <- data.frame(x = rnorm(40))
dd$y <- rnorm(40, 1 + 0.5 * dd$x)
fit <- q(frm(bf(y ~ x), family = gaussian(), data = dd))
tpl <- fit$frame$par_template
est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
M <- matrix(rep(est, each = 300) + rnorm(300 * length(est), 0, 0.05), 300,
            dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
ds <- structure(list(stanfit = NULL, draws = M, fit = fit),
                class = "frmtmb_draws")
tc <- function(e) tryCatch({
  v <- q(e); paste("ok, dim", paste(dim(v), collapse = " x "))
}, error = function(err) paste("ERROR:", substr(conditionMessage(err), 1, 120)))
cat("psis(ds):                ", tc(psis(ds)), "\n")
cat("psis(ds, ndraws = 100):  ", tc(psis(ds, ndraws = 100)), "\n")
cat("psis(ds, draw_ids = 1:100):", tc(psis(ds, draw_ids = 1:100)), "\n")
cat("psis(ds, junk = 1):      ", tc(psis(ds, junk = 1)), "\n")
