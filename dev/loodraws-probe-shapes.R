# Shapes of the draws methods the loo_*() functions build on, on the
# models the brief names. Log: dev/loodraws-log/probe-shapes.txt
.libPaths(c("C:/Users/adf44/source/r/wt-loodraws-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-loodraws/dev/stan-cache")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", dirname(system.file(package = "frmtmb.sample")), "\n")
set.seed(1)
n <- 40
dd <- data.frame(x = rnorm(n), g = factor(rep(1:8, 5)))
dd$y <- rnorm(n, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g])
dd$y2 <- rnorm(n, -1 + 0.3 * dd$x)
dd$k <- rpois(n, exp(0.2 + 0.4 * dd$x))
dd$o <- factor(sample(1:4, n, TRUE), ordered = TRUE)
smp <- function(...) suppressWarnings(suppressMessages(
  frm_sample(..., chains = 2, iter = 200, warmup = 100, seed = 1,
             refresh = 0)))
t0 <- proc.time()
ds_mv <- smp(mvbf(bf(y ~ x), bf(y2 ~ x), rescor = FALSE),
             family = gaussian(), data = dd)
ds_or <- smp(bf(o ~ x), family = sratio(), data = dd)
print(proc.time() - t0)
show <- function(lab, expr) {
  v <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(lab, ": ", if (is.character(v)) v else paste(dim(v), collapse = " x "),
      "\n", sep = "")
}
show("mv ll", log_lik(ds_mv))
show("mv ll resp y", log_lik(ds_mv, resp = "y"))
show("mv epred", posterior_epred(ds_mv))
show("mv linpred", posterior_linpred(ds_mv))
show("mv predict", posterior_predict(ds_mv))
show("mv epred y", posterior_epred(ds_mv, resp = "y"))
show("or ll", log_lik(ds_or))
show("or epred", posterior_epred(ds_or))
show("or linpred", posterior_linpred(ds_or))
show("or linpred thres", posterior_linpred(ds_or, incl_thres = TRUE))
show("or predict", posterior_predict(ds_or))
show("linpred scale", posterior_linpred(ds_or, scale = "response"))
show("linpred junk", posterior_linpred(ds_or, junk = 1))
show("epred junk", posterior_epred(ds_or, junk = 1))
print(dimnames(posterior_epred(ds_mv))[[3]])
print(attributes(log_lik(ds_mv))[c("chain_id")] |> lengths())
