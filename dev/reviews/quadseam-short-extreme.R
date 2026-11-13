# Short check of lane quadseam's fix round, item 2: fixed-effects fits
# (no random effects, so the fit-end headroom reads obj$he()) that stop
# where the AD Hessian may not be finite: a poisson all-zero response,
# complete separation, all rows right-censored, a huge-count poisson, and
# a binomial at p = 1. Per case: is he() finite at the stop, the
# warnings and messages by class and text, the convergence code, and
# seconds.
#   Rscript quadseam-short-extreme.R <base|lane>
a <- commandArgs(TRUE)
arm <- a[1]
base <- "C:/Users/adf44/source/r/rellib-r20"
lane <- "C:/Users/adf44/source/r/wt-quadseam-lib"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (arm == "lane") c(lane, base, ul) else c(base, ul))
suppressMessages(library(frmtmb))
cat("lib:", arm, dirname(find.package("frmtmb")), "\n")
set.seed(7)
n <- 60
x <- rnorm(n)
cases <- list(
  pois_zero = list(d = data.frame(y = 0L, x = x), f = bf(y ~ x),
                   fam = poisson()),
  sep = list(d = data.frame(y = as.integer(x > 0), x = x), f = bf(y ~ x),
             fam = bernoulli()),
  sep_quasi = list(d = data.frame(y = as.integer(x > 0.3), x = x,
                                  z = rnorm(n)),
                   f = bf(y ~ x + z), fam = bernoulli()),
  cens_all = list(d = data.frame(y = 1 + x, cc = "right", x = x),
                  f = bf(y | cens(cc) ~ x), fam = gaussian()),
  pois_big = list(d = data.frame(y = rpois(n, exp(12 + 0.1 * x)), x = x),
                  f = bf(y ~ x), fam = poisson()),
  bin_one = list(d = data.frame(k = 10L, N = 10L, x = x),
                 f = bf(k | trials(N) ~ x), fam = binomial()),
  gamma_ok = list(d = data.frame(y = rgamma(n, 3, 3 / exp(1 + 0.4 * x)),
                                 x = x), f = bf(y ~ x),
                  fam = Gamma(link = "log"))
)
for (nm in names(cases)) {
  cs <- cases[[nm]]
  w <- character(); ms <- character()
  t0 <- proc.time()[[3]]
  f <- tryCatch(withCallingHandlers(
    frm(cs$f, family = cs$fam, data = cs$d),
    warning = function(x) {
      w <<- c(w, paste0(class(x)[1], ": ",
                        substr(gsub("\\s+", " ", conditionMessage(x)), 1, 110)))
      invokeRestart("muffleWarning")
    },
    message = function(x) {
      ms <<- c(ms, class(x)[1]); invokeRestart("muffleMessage")
    }), error = function(e) e)
  tt <- proc.time()[[3]] - t0
  if (inherits(f, "error")) {
    cat(sprintf("CASE %-9s %s ERROR %s\n", nm, arm,
                substr(conditionMessage(f), 1, 150)))
    next
  }
  H <- tryCatch(f$obj$he(f$opt$par), error = function(e) NULL)
  cat(sprintf("CASE %-9s %s he_finite %s conv %d ll %.6f secs %.2f msg %s\n",
              nm, arm, if (is.null(H)) "NA" else all(is.finite(H)),
              f$opt$convergence, as.numeric(logLik(f)), tt,
              if (length(ms)) paste(ms, collapse = ",") else "-"))
  for (x in w) cat("    W", x, "\n")
}
cat("END\n")
