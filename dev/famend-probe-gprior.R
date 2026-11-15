# Lane famend: 3560 probe. A positive-support prior on a coefficient with
# no natural lower bound: what the start, the objective and lb = 0 do.
#   Rscript dev/famend-probe-gprior.R <lib or "base">
a <- commandArgs(trailingOnly = TRUE)
base <- "C:/Users/adf44/source/r/rellib-r21"
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (a[1] != "base") a[1], base, ulib))
suppressPackageStartupMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "\n")
cap <- function(expr) {
  w <- character()
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }, message = function(x) invokeRestart("muffleMessage")),
    error = function(e) structure(conditionMessage(e), class = "err"))
  list(v = v, w = w)
}
show <- function(lab, r) {
  cat("\n==", lab, "\n")
  for (x in r$w) cat("   W:", substr(gsub("\\s+", " ", x), 1, 160), "\n")
  if (inherits(r$v, "err")) {
    cat("   ERR:", substr(r$v, 1, 300), "\n")
    return(invisible())
  }
  f <- r$v
  cat("   code", f$opt$convergence, "objective", f$opt$objective,
      "betad", format(f$estimates$betad, digits = 6), "\n")
}
set.seed(20261021)
n <- 80
d <- data.frame(x = rnorm(n), fac = factor(rep(c("a", "b"), 40)))
eta <- 0.5 + 0.4 * d$x + rnorm(n, 0, 0.5)
d$yg <- rgamma(n, shape = 3, rate = 3 / exp(eta))
d$yn <- rgamma(n, shape = 0.4, rate = 0.4 / exp(eta))
fo <- bf(yg ~ x, shape ~ 0 + fac)
P <- function(...) set_prior(...)
st <- par_template(fo, data = d, family = Gamma("log"))
cat("start betad:", format(st$betad), "\n")
show("gamma(2, 1), no lb", cap(frm(fo, family = Gamma("log"), data = d,
                                   prior = P("gamma(2, 1)", class = "b",
                                             dpar = "shape"))))
show("gamma(2, 1), lb = 0", cap(frm(fo, family = Gamma("log"), data = d,
                                    prior = P("gamma(2, 1)", class = "b",
                                              dpar = "shape", lb = 0))))
show("gamma(2, 1), no lb, start betad = 1", cap(frm(fo, family = Gamma("log"),
  data = d, prior = P("gamma(2, 1)", class = "b", dpar = "shape"),
  start = list(betad = c(1, 1)))))
show("lognormal(0, 1), no lb", cap(frm(fo, family = Gamma("log"), data = d,
  prior = P("lognormal(0, 1)", class = "b", dpar = "shape"))))
show("exponential(1), no lb", cap(frm(fo, family = Gamma("log"), data = d,
  prior = P("exponential(1)", class = "b", dpar = "shape"))))
show("normal(0, 1) (absent)", cap(frm(fo, family = Gamma("log"), data = d,
  prior = P("normal(0, 1)", class = "b", dpar = "shape"))))
fn <- bf(yn ~ x, shape ~ 0 + fac)
show("shape 0.4 data: gamma(2, 1), lb = 0", cap(frm(fn, family = Gamma("log"),
  data = d, prior = P("gamma(2, 1)", class = "b", dpar = "shape", lb = 0))))
show("shape 0.4 data: gamma(2, 1), start 1", cap(frm(fn,
  family = Gamma("log"), data = d, start = list(betad = c(1, 1)),
  prior = P("gamma(2, 1)", class = "b", dpar = "shape"))))
# on a mu coefficient
show("gaussian: gamma(2, 1) on b of mu, no lb", cap(frm(yg ~ x, data = d,
  prior = P("gamma(2, 1)", class = "b"))))
show("gaussian: gamma(2, 1) on b of mu, lb = 0", cap(frm(yg ~ x, data = d,
  prior = P("gamma(2, 1)", class = "b", lb = 0))))
show("gaussian: uniform(0, 3) on b of mu, no bounds", cap(frm(yg ~ x,
  data = d, prior = P("uniform(0, 3)", class = "b"))))
show("gaussian: beta(2, 2) on b of mu, no bounds", cap(frm(yg ~ x,
  data = d, prior = P("beta(2, 2)", class = "b"))))
cat("\ndone\n")
