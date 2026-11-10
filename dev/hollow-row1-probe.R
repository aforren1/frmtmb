# Lane hollow, Row 1: what frmtmb's hypothesis() reports for the
# evidence ratio on (a) an ML fit, (b) a MAP fit with proper priors and
# (c) draws with prior draws, and what a Laplace (normal) approximation
# of the posterior at the MAP would give in place of (b)'s NA.
#   Rscript dev/hollow-row1-probe.R [lib] > dev/hollow-out/row1-probe.txt
# lib defaults to the base build, rellib-r19.
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-hollow/dev/stan-cache")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
cat("lib:", find.package("frmtmb.sample"),
    format(packageVersion("frmtmb.sample")), "\n")
q <- function(expr) suppressWarnings(suppressMessages(expr))

set.seed(394)
n <- 120
dd <- data.frame(x = rnorm(n), g = factor(rep(1:12, 10)))
dd$y <- rnorm(n, 1 + 0.15 * dd$x + rnorm(12, 0, 0.6)[dd$g], 1)
pri <- set_prior("normal(0, 1)", class = "b") +
  set_prior("cauchy(0, 2)", class = "sd", lb = 0)
hs <- c("x = 0", "x > 0")

cat("\n== (a) ML fit, flat priors ==\n")
ml <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
ha <- hypothesis(ml, hs)
print(ha$hypothesis, digits = 6)
hsd <- hypothesis(ml, "Intercept = 0", class = "sd", group = "g")
print(hsd$hypothesis, digits = 6)
print(hsd)

cat("\n== (b) MAP fit, normal(0, 1) on b, cauchy(0, 2) on sd ==\n")
map <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
           prior = pri)
hb <- hypothesis(map, hs)
print(hb$hypothesis, digits = 6)
hbsd <- hypothesis(map, "Intercept = 0", class = "sd", group = "g")
print(hbsd$hypothesis, digits = 6)
cat("-- print() of the MAP hypothesis --\n")
print(hbsd)

cat("\n-- Laplace (normal at the MAP) stand-ins for (b)'s NA --\n")
ex <- hb$hypothesis$Estimate[1]; se <- hb$hypothesis$Est.Error[1]
lap_x <- dnorm(0, ex, se) / dnorm(0, 0, 1)
cat("x = 0, normal on b_x's own scale:", format(lap_x, digits = 6), "\n")
cat("x > 0, normal posterior odds:",
    format(pnorm(ex / se) / pnorm(-ex / se), digits = 6), "\n")
es <- hbsd$hypothesis$Estimate[1]; ss <- hbsd$hypothesis$Est.Error[1]
half_cauchy0 <- 2 * dcauchy(0, 0, 2)
cat("sd = 0, normal on the sd scale:",
    format(dnorm(0, es, ss) / half_cauchy0, digits = 6), "\n")
cat("sd = 0, normal on the log-sd scale (the optimizer's): density of ",
    "sd at 0 is 0, ratio 0\n", sep = "")

cat("\n== (c) draws of the MAP fit, sample_prior = 'yes' ==\n")
ds <- q(frm_sample(map, chains = 4, iter = 2000, warmup = 1000,
                   seed = 3940, refresh = 0, sample_prior = "yes"))
hc <- hypothesis(ds, hs)
print(hc$hypothesis, digits = 6)
hcsd <- hypothesis(ds, "Intercept = 0", class = "sd", group = "g")
print(hcsd$hypothesis, digits = 6)
cat("prior_draws columns:", paste(names(prior_draws(ds)), collapse = " "),
    "\n")
pd <- hcsd$prior_samples$H1
cat("prior draws of sd_g__Intercept: n", length(pd), " min",
    format(min(pd), digits = 4), " KS p vs half-Cauchy(0, 2):",
    format(ks.test(pd, function(v) 2 * pcauchy(v, 0, 2) - 1)$p.value,
           digits = 4), "\n")
if (requireNamespace("brms", quietly = TRUE)) {
  dr <- utils::getFromNamespace("density_ratio", "brms")
  for (k in list(list("x = 0", hc, 1L), list("sd = 0", hcsd, 1L))) {
    h <- k[[2]]
    b <- dr(h$samples[[k[[3]]]], h$prior_samples[[k[[3]]]], point = 0)
    cat(k[[1]], ": frmtmb Evid.Ratio", format(h$hypothesis$Evid.Ratio[1],
                                              digits = 10),
        " brms:::density_ratio() on the same draws", format(b, digits = 10),
        " identical:", identical(h$hypothesis$Evid.Ratio[1], b), "\n")
  }
}
cat("ROW1 PROBE END\n")
