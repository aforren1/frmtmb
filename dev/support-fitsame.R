# Lane support, part A: the fit is unchanged. Each model below is fitted
# and its estimates and logLik() printed in hexadecimal, so that two
# builds can be compared bit for bit; run it once per build, pinned to
# one core (dev/rtmb-pitfalls.md, item 21), and diff the outputs.
# Then the per-row log-likelihood a reader of draws computes,
# frmtmb.sample's draws_row_loglik() on doubles, is summed at the
# estimates and set against logLik(), the tape's value.
#
#   Rscript dev/support-fitsame.R   (FRMTMB_PORT_LIB picks the build)
lib <- Sys.getenv("FRMTMB_PORT_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", dirname(system.file(package = "frmtmb")), "frmtmb",
    format(packageVersion("frmtmb")), "\n")
rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
set.seed(31)
n <- 600
d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
nu <- 1 + exp(log(3) + 0.2 * d$z)
d$yf <- exp(1 + 0.3 * d$x) / gamma(1 - 1 / nu) * stats::rexp(n)^(-1 / nu)
d$yg <- 1 + 0.5 * d$x + exp(-0.3 + 0.2 * d$z) *
  (stats::rexp(n)^(-0.15) - 1) / 0.15
d$ye <- stats::rexp(n, 1 / exp(0.5 + 0.3 * d$x))
# two GEV components with xi 0.3, the second's support starting at
# 6 - 0.5 / 0.3 = 4.33: most rows of the first lie below it, so at the
# fit that component misses them and adds exp(-Inf) = 0 there
set.seed(32)
m <- 400
dm <- data.frame(y = c(rg(m / 2, 0, 1, 0.3), rg(m / 2, 6, 0.5, 0.3)))
# censoring and truncation: right-censor the top decile, left-censor the
# bottom one, truncate below at a point under every row
cen <- function(y) ifelse(y > stats::quantile(y, 0.9), 1,
                          ifelse(y < stats::quantile(y, 0.1), -1, 0))
d$cg <- cen(d$yg)
d$cf <- cen(d$yf)
d$ce <- cen(d$ye)
d$lbg <- min(d$yg) - 1
d$lbf <- min(d$yf) / 2
models <- list(
  gev = list(bf(yg ~ x, sigma ~ z), gen_extreme_value(), d),
  gev_cens_trunc = list(bf(yg | cens(cg) + trunc(lb = lbg) ~ x),
                        gen_extreme_value(), d),
  gev_mixture = list(bf(y ~ 1), mixture(gen_extreme_value(),
                                        gen_extreme_value()), dm),
  frechet = list(bf(yf ~ x, nu ~ z), frechet(), d),
  frechet_cens_trunc = list(bf(yf | cens(cf) + trunc(lb = lbf) ~ x),
                            frechet(), d),
  weibull_cens = list(bf(ye | cens(ce) ~ x), weibull(), d),
  exponential_cens = list(bf(ye | cens(ce) ~ x), exponential(), d)
)
hex <- function(v) paste(sprintf("%a", v), collapse = " ")
for (nm in names(models)) {
  mm <- models[[nm]]
  fit <- tryCatch(withCallingHandlers(
    frm(mm[[1]], family = mm[[2]], data = mm[[3]]),
    warning = function(w) {
      cat("  warning:", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    }), error = function(e) conditionMessage(e))
  if (is.character(fit)) {
    cat(nm, "ERROR", fit, "\n")
    next
  }
  ll <- as.numeric(logLik(fit))
  cat(nm, "par", hex(fit$opt$par), "\n")
  cat(nm, "logLik", hex(ll), sprintf("(%.12g)", ll), "\n")
  rows <- tryCatch(frmtmb.sample:::draws_row_loglik(fit, NULL),
                   error = function(e) conditionMessage(e))
  if (is.character(rows)) {
    cat(nm, "reader ERROR", rows, "\n")
  } else {
    cat(sprintf("%s reader: %d rows, %d not finite, sum - logLik = %.3g\n",
                nm, length(rows), sum(!is.finite(rows)), sum(rows) - ll))
  }
}
