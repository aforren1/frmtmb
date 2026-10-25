# The cox() fits that do not use the new arguments are unchanged:
# logLik, estimates, vcov and the baseline weights of each fit, written
# to an rds per arm and compared with identical() by
# dev/bhaz-bitwise-cmp.R.
#   Rscript dev/bhaz-bitwise.R <lane|base> <out.rds>
args <- commandArgs(trailingOnly = TRUE)
arm <- args[[1L]]
libs <- c("C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-bhaz-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb.survival))
cat("arm", arm, "frmtmb at", find.package("frmtmb"), "survival at",
    find.package("frmtmb.survival"), "\n")
set.seed(33)
n <- 400
x <- rnorm(n)
g <- factor(sample(1:20, n, TRUE))
u <- rnorm(20, 0, 0.5)[g]
tt <- stats::rweibull(n, 1.4, exp(-(-0.5 + 0.8 * x + u) / 1.4))
ct <- stats::rexp(n, 0.3)
dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x,
                 g = g, entry = pmin(tt, ct) * runif(n, 0, 0.3))
dd$upper <- ifelse(dd$cens == 1, dd$time * 1.5, dd$time)
dd$c2 <- ifelse(dd$cens == 1, 2, 0)
fits <- list(
  default = function() frm(bf(time | cens(cens) ~ x), family = cox(),
                           data = dd),
  df6 = function() frm(bf(time | cens(cens) ~ x), family = cox(df = 6),
                       data = dd),
  df9 = function() frm(bf(time | cens(cens) ~ x), family = cox(df = 9),
                       data = dd),
  deg2 = function() frm(bf(time | cens(cens) ~ x),
                        family = cox(df = 5, degree = 2, intercept = FALSE),
                        data = dd),
  identity = function() frm(bf(time | cens(cens) ~ x),
                            family = cox(link = "identity"),
                            data = transform(dd, x = abs(x))),
  frailty = function() frm(bf(time | cens(cens) ~ x + (1 | g)),
                           family = cox(), data = dd),
  trunc = function() frm(bf(time | cens(cens) + trunc(lb = entry) ~ x),
                         family = cox(), data = dd),
  interval = function() frm(bf(time | cens(c2, upper) ~ x), family = cox(),
                            data = dd),
  string = function() frm(bf(time | cens(cens) ~ x), family = "cox",
                          data = dd)
)
out <- lapply(names(fits), function(nm) {
  f <- tryCatch(suppressWarnings(fits[[nm]]()), error = function(e) e)
  if (inherits(f, "error")) {
    cat(nm, "ERROR", conditionMessage(f), "\n")
    return(list(error = conditionMessage(f)))
  }
  cat(sprintf("%-9s logLik %.15g\n", nm, as.numeric(logLik(f))))
  list(logLik = logLik(f), estimates = f$estimates, vcov = vcov(f),
       baseline = cox_baseline(f))
})
names(out) <- names(fits)
saveRDS(out, args[[2L]])
