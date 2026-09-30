# fn and gr in hexadecimal at a fixed perturbed point, for models this
# lane's changes pass through without adding a term: the count families
# whose densities were rewritten for rate(), the frame's per-response
# loops (a multivariate model, mi(), cs(), a threshold-only ordinal
# model), and predictions. Run once per library and diff:
#   Rscript dev/aterms2-bitwise.R lane > dev/aterms2-log-bitwise-lane.txt
#   Rscript dev/aterms2-bitwise.R base > dev/aterms2-log-bitwise-base.txt
arm <- commandArgs(TRUE)[[1]]
libs <- c("C:/Users/adf44/source/r/wt-aterms2-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
hex <- function(x) sprintf("%a", as.numeric(x))
set.seed(31)
n <- 150
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, 15)),
                w = runif(n, 0.5, 2))
d$yc <- rpois(n, exp(0.3 + 0.4 * d$x + rnorm(10, 0, 0.5)[d$g]))
d$yn <- rnbinom(n, mu = exp(0.3 + 0.4 * d$x), size = 2)
d$y1 <- 1 + d$x + rnorm(n)
d$y2 <- 0.5 - d$z + rnorm(n)
d$o <- factor(sample(1:4, n, TRUE), ordered = TRUE)
d$ym <- d$y1
d$ym[c(3, 9, 40)] <- NA
d$cc <- rep(c(0, 0, 1), n / 3)
models <- list(
  poisson_re = function() frm(yc ~ x + (1 | g), data = d, family = poisson()),
  poisson_ident = function() frm(yc ~ x, data = d,
                                 family = poisson("sqrt")),
  poisson_cens = function() frm(yc | cens(cc) ~ x, data = d,
                                family = poisson()),
  negbin_shape = function() frm(bf(yn ~ x, shape ~ z), data = d,
                                family = negbinomial()),
  negbin_ident = function() frm(yn ~ 1, data = d,
                                family = negbinomial("identity")),
  geometric = function() frm(yn ~ x, data = d, family = geometric()),
  mv_re = function() frm(bf(y1 ~ x + (1 | p | g)) + bf(y2 ~ z + (1 | p | g)),
                         data = d, family = gaussian()),
  mi = function() frm(bf(y2 ~ mi(ym) + z) + bf(ym | mi() ~ x), data = d,
                      family = gaussian()),
  ordinal_cs = function() frm(o ~ cs(x), data = d, family = sratio()),
  ordinal_only = function() frm(o ~ 1 + (1 | g), data = d,
                                family = cumulative())
)
for (nm in names(models)) {
  f <- suppressWarnings(suppressMessages(models[[nm]]()))
  set.seed(7)
  p <- f$opt$par + rnorm(length(f$opt$par), 0, 0.05)
  cat(nm, "fn", hex(f$obj$fn(p)), "\n")
  cat(nm, "gr", hex(f$obj$gr(p)), "\n")
  cat(nm, "logLik", hex(logLik(f)), "\n")
  if (length(f$spec$responses) == 1L) {
    fv <- suppressWarnings(fitted(f))
    cat(nm, "fitted", hex(head(as.numeric(fv), 20)), "\n")
    if (!identical(f$spec$responses[[1]]$family$type, "ordinal")) {
      cat(nm, "pearson", hex(head(residuals(f, type = "pearson")[, 1], 10)),
          "\n")
    }
  }
}
