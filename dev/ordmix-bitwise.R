# Plain ordinal and mixture fits the lane must leave unchanged to the
# bit: opt$par, logLik, fitted(), simulate(seed = 1), variables(),
# fixef() and print(). Data seed 20261050.
# Usage: Rscript dev/ordmix-bitwise.R <base|lane>  (writes
# dev/ordmix-bitwise-<arm>.rds), then Rscript dev/ordmix-bitwise.R
# compare.
args <- commandArgs(TRUE)
arm <- args[1]
root <- "C:/Users/adf44/source/r/frmtmb-wt-ordmix/dev"
if (arm == "compare") {
  a <- readRDS(file.path(root, "ordmix-bitwise-base.rds"))
  b <- readRDS(file.path(root, "ordmix-bitwise-lane.rds"))
  n <- 0L
  same <- 0L
  for (m in names(a)) {
    for (o in names(a[[m]])) {
      n <- n + 1L
      ok <- identical(a[[m]][[o]], b[[m]][[o]])
      same <- same + ok
      if (!ok) cat("DIFFERS", m, o, "\n")
    }
  }
  cat(sprintf("BITWISE %d of %d outputs identical\n", same, n))
  quit(save = "no")
}
libs <- c("C:/Users/adf44/source/r/rellib-r5",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ordmix-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, find.package("frmtmb"), "\n")
set.seed(20261050)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
u <- rlogis(n) / exp(0.3 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -1) + (u > 0.3) + (u > 1.5)
d$yh <- ifelse(runif(n) < 0.2, 0L, d$y)
d$yn <- rnorm(n, ifelse(runif(n) < 0.4, 2, -1) + 0.5 * d$x)
models <- list(
  cum = function() frm(bf(y ~ x), family = cumulative(), data = d),
  cum_probit_gr = function() frm(bf(y | thres(gr = g) ~ x),
                                 family = cumulative("probit"), data = d),
  cum_disc = function() frm(bf(y ~ x, disc ~ 0 + z), family = cumulative(),
                            data = d),
  sratio_cs = function() frm(bf(y ~ cs(x)), family = sratio(), data = d),
  acat_gr = function() frm(bf(y | thres(gr = g) ~ x), family = acat(),
                           data = d),
  cratio_stz = function() frm(bf(y ~ x),
                              family = cratio(threshold = "sum_to_zero"),
                              data = d),
  hurdle = function() frm(bf(yh ~ x, hu ~ z), family = hurdle_cumulative(),
                          data = d),
  hurdle_equi = function() frm(bf(yh ~ x),
                               family = hurdle_cumulative(
                                 threshold = "equidistant"), data = d),
  gauss_mix = function() frm(bf(yn ~ x), family = mixture(gaussian(),
                                                          gaussian()),
                             data = d)
)
out <- list()
for (m in names(models)) {
  fit <- suppressWarnings(models[[m]]())
  out[[m]] <- list(
    par = fit$opt$par, logLik = as.numeric(logLik(fit)),
    fitted = fitted(fit), sim = simulate(fit, nsim = 2, seed = 1),
    variables = variables(fit), fixef = fixef(fit),
    print = capture.output(print(fit)))
}
saveRDS(out, file.path(root, paste0("ordmix-bitwise-", arm, ".rds")))
