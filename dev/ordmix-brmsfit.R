# Fits against brms 2.23.0 with seeds: brms's posterior (its default
# priors, one chain started at frmtmb's estimates so that both describe
# the same labeling of the components) against frmtmb's maximum
# likelihood fit, row by row of fixef(), in units of brms's posterior
# standard deviation; and brms's fixef() row names and order against
# frmtmb's. Usage: Rscript dev/ordmix-brmsfit.R <case>
# Output: dev/ordmix-brmsfit-log/<case>.txt (driver ordmix-brmsfit.ps1)
args <- commandArgs(TRUE)
case <- if (length(args)) args[1] else "none"
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(brms)
})
cat("case", case, "frmtmb", find.package("frmtmb"), "brms",
    format(packageVersion("brms")), "\n")
seed <- 20261040 + match(case, c("none", "mu", "hgr", "hcs", "none2"))
set.seed(seed)
n <- 800
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
cls <- rbinom(n, 1, 0.4)
lat <- switch(case,
  mu = ifelse(cls == 1, 2 * d$x, -1.5 * d$x),
  ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1)) + rlogis(n)
d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
u <- rlogis(n) + 0.8 * d$x
d$yh <- ifelse(runif(n) < plogis(-1 + 0.5 * d$z), 0L,
               1L + (u > -1 + 0.1 * d$x) + (u > ifelse(d$g == "a", 0.3, 0.6)) +
                 (u > 1.5 - 0.1 * d$x))
spec <- switch(case,
  none = list(ff = frmtmb::bf(y ~ x), fb = brms::bf(y ~ x),
              famf = frmtmb::mixture(frmtmb::cumulative(), frmtmb::sratio()),
              famb = brms::mixture(brms::cumulative(), brms::sratio())),
  none2 = list(ff = frmtmb::bf(y ~ x), fb = brms::bf(y ~ x),
               famf = frmtmb::mixture(frmtmb::cumulative(),
                                      frmtmb::cumulative()),
               famb = brms::mixture(brms::cumulative(), brms::cumulative())),
  mu = list(ff = frmtmb::bf(y ~ x), fb = brms::bf(y ~ x),
            famf = frmtmb::mixture(frmtmb::cumulative(), frmtmb::cumulative(),
                                   order = "mu"),
            famb = brms::mixture(brms::cumulative(), brms::cumulative(),
                                 order = "mu")),
  hgr = list(ff = frmtmb::bf(yh | thres(gr = g) ~ x, hu ~ z),
             fb = brms::bf(yh | thres(gr = g) ~ x, hu ~ z),
             famf = frmtmb::hurdle_cumulative(),
             famb = brms::hurdle_cumulative()),
  hcs = list(ff = frmtmb::bf(yh ~ cs(x)), fb = brms::bf(yh ~ cs(x)),
             famf = frmtmb::hurdle_cumulative("probit"),
             famb = brms::hurdle_cumulative("probit")))
fit <- frm(spec$ff, family = spec$famf, data = d)
fx <- fixef(fit)
print(fx)
# brms's parameters at frmtmb's estimates, for the chain's start
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-ordmix/tests/testthat/test-ordinal-mixture-brms.R")
a <- grep("^omb_par_names <- function", src)
b <- grep("^omb_check <- function", src) - 1L
eval(parse(text = src[a:b]))
environment(omb_translate) <- asNamespace("frmtmb")
code <- suppressWarnings(as.character(
  brms::stancode(spec$fb, data = d, family = spec$famb)))
sdat <- suppressWarnings(brms::standata(spec$fb, data = d,
                                        family = spec$famb))
init <- omb_translate(fit, fit$estimates, code, sdat)
bf_fit <- suppressWarnings(brms::brm(spec$fb, data = d, family = spec$famb,
                                     chains = 1, iter = 1500, warmup = 750,
                                     seed = 1, refresh = 0,
                                     init = list(init),
                                     backend = "rstan"))
bx <- brms::fixef(bf_fit)
print(bx)
cat("ROWS identical names and order:", identical(rownames(bx), rownames(fx)),
    "\n")
cat("ROWS brms:", rownames(bx), "\n")
cat("ROWS frmtmb:", rownames(fx), "\n")
common <- intersect(rownames(bx), rownames(fx))
z <- (fx[common, "Estimate"] - bx[common, "Estimate"]) / bx[common, "Est.Error"]
se_ratio <- fx[common, "Est.Error"] / bx[common, "Est.Error"]
for (r in common) {
  cat(sprintf("ROW %s case=%s frmtmb=%.4f brms=%.4f z=%.3f se_ratio=%.3f\n",
              r, case, fx[r, "Estimate"], bx[r, "Estimate"], z[[r]],
              se_ratio[[r]]))
}
cat(sprintf("SUMMARY case=%s rows=%d max|z|=%.3f se_ratio_range=[%.3f, %.3f] rhat_max=%.3f\n",
            case, length(common), max(abs(z)), min(se_ratio),
            max(se_ratio), max(brms::rhat(bf_fit), na.rm = TRUE)))
