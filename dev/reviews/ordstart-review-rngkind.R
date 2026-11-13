# Review of lane ordstart: the private stream is set.seed(20261111) under
# the session's RNG kind. Does a session with another kind (parallel
# work sets L'Ecuyer-CMRG; R < 3.6 scripts set sample.kind = "Rounding",
# which does not touch rnorm) get another fit from the same call? 40-fit
# data, all 40 seeds.
#   Rscript ordstart-review-rngkind.R
.libPaths(c("C:/Users/adf44/source/r/wt-ordstart-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
q <- function(e) suppressWarnings(suppressMessages(e))
z1 <- frmtmb:::mixture_starts_draws(3, 1)
RNGkind("L'Ecuyer-CMRG")
z2 <- frmtmb:::mixture_starts_draws(3, 1)
RNGkind("default", "default", "default")
cat("draws, default kind:", round(z1, 4), "| L'Ecuyer-CMRG:", round(z2, 4),
    "\n")
nd <- 0L
for (nm in c("cum_cum", "cum_sratio")) for (s in 1:20) {
  fam <- switch(nm, cum_cum = mixture(cumulative(), cumulative()),
                cum_sratio = mixture(cumulative(), sratio()))
  d <- mk(s, mix = TRUE)
  RNGkind("default", "default", "default")
  a <- as.numeric(logLik(q(frm(bf(y ~ x + z), family = fam, data = d))))
  RNGkind("L'Ecuyer-CMRG")
  b <- as.numeric(logLik(q(frm(bf(y ~ x + z), family = fam, data = d))))
  RNGkind("default", "default", "default")
  if (abs(a - b) > 1e-6) nd <- nd + 1L
  cat(sprintf("RNGKIND %s %2d default %.4f lecuyer %.4f diff %.4f\n", nm, s,
              a, b, b - a))
}
cat("fits whose logLik differs by RNG kind:", nd, "of 40\n")
