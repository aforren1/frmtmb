# Lane ordstarts, build round: every fit without an ordinal mixture is
# bitwise unchanged. Fits a set of models under one arm (lane: the lane
# library first; base: rellib-r22) and saves each fit's optimum, logLik,
# code and estimates to an RDS; compare the two with
#   Rscript dev/ordstarts-bitwise.R <lane|base> <out.rds>
#   Rscript dev/ordstarts-bitwise.R compare <lane.rds> <base.rds>
# The models: gaussian, poisson with (1 | g), cumulative, sratio with
# cs(), a gaussian mixture, a negbinomial hurdle, a beta, and two
# ordinal mixtures with starts = 0 (the path the swap and the default
# starts leave as it was).
a <- commandArgs(trailingOnly = TRUE)
if (a[1] == "compare") {
  L <- readRDS(a[2])
  B <- readRDS(a[3])
  stopifnot(identical(names(L), names(B)))
  for (nm in names(L)) {
    same <- identical(L[[nm]], B[[nm]])
    cat(sprintf("%-28s identical %s\n", nm, same))
  }
  cat("ALL IDENTICAL", all(mapply(identical, L, B)), "\n")
  quit(save = "no")
}
arm <- a[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-ordstarts-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
set.seed(42)
n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep_len(1:20, n)))
u <- rnorm(20, 0, 0.5)
d$yg <- 1 + 0.5 * d$x + u[d$g] + rnorm(n)
d$yp <- rpois(n, exp(0.3 + 0.4 * d$x + u[d$g]))
lat <- 0.8 * d$x + rlogis(n)
d$yo <- 1L + (lat > -1) + (lat > 0) + (lat > 1)
d$ym <- ifelse(rbinom(n, 1, 0.4) == 1, rnorm(n, 2 + d$x), rnorm(n, -1))
d$yh <- ifelse(rbinom(n, 1, 0.3) == 1, 0,
               pmax(1, rnbinom(n, mu = exp(1 + 0.3 * d$x), size = 2)))
d$yb <- plogis(0.3 * d$x + rnorm(n, 0, 0.5))
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
dm <- mk(13, n = 300, mix = TRUE)
q <- function(e) suppressWarnings(suppressMessages(e))
keep <- function(f) {
  list(par = f$opt$par, obj = f$opt$objective, code = f$opt$convergence,
       ll = as.numeric(logLik(f)), est = f$estimates[setdiff(
         names(f$estimates), character(0))])
}
fits <- list(
  gaussian = function() frm(bf(yg ~ x + (1 | g)), data = d),
  poisson_re = function() frm(bf(yp ~ x + (1 | g)), family = poisson(),
                              data = d),
  cumulative = function() frm(bf(yo ~ x + z), family = cumulative(),
                              data = d),
  sratio_cs = function() frm(bf(yo ~ cs(x) + z), family = sratio(),
                             data = d),
  mix_gauss = function() frm(bf(ym ~ x), family = mixture(gaussian(),
                                                         gaussian()),
                             data = d),
  hurdle_nb = function() frm(bf(yh ~ x), family = hurdle_negbinomial(),
                             data = d),
  beta = function() frm(bf(yb ~ x), family = Beta(), data = d),
  ordmix_cc_starts0 = function() {
    frm(bf(y ~ x + z), family = mixture(cumulative(), cumulative()),
        data = dm, control = frmtmb_control(starts = 0))
  },
  ordmix_cs_starts0 = function() {
    frm(bf(y ~ x + z), family = mixture(cumulative(), sratio()), data = dm,
        control = frmtmb_control(starts = 0))
  })
res <- lapply(fits, function(fn) {
  f <- tryCatch(q(fn()), error = function(e) conditionMessage(e))
  if (is.character(f)) f else keep(f)
})
saveRDS(res, a[2])
cat("saved", length(res), "fits\n")
