# Lane optima, punch round 1, B1: a model without mo() samples as it
# did on 0.68.1. Three models (gaussian with (1 | g), poisson, a
# cumulative with cs()), frm_sample(chains = 2, iter = 600, seed = 1)
# with the default priors and with prior = "flat"; the draws matrices
# are saved per arm and compared by dev/optima-p1-nomo-cmp.R. Run each
# arm pinned to one core (cmd start /affinity), dev/rtmb-pitfalls.md
# item 21.
#   Rscript dev/optima-p1-nomo.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
setwd("C:/Users/adf44/source/r/frmtmb-wt-optima")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(3)
g <- factor(rep(1:12, each = 10))
d <- data.frame(g = g, x = rnorm(120), z = rnorm(120))
d$y <- 1 + 0.5 * d$x + rnorm(12, 0, 0.7)[g] + rnorm(120)
d$c <- rpois(120, exp(0.3 + 0.4 * d$x))
d$o <- as.integer(cut(d$x + rnorm(120), c(-Inf, -1, 0, 1, Inf)))
fits <- list(
  gauss_re = frm(y ~ x + (1 | g), data = d),
  pois = frm(c ~ x, family = poisson(), data = d),
  cum_cs = suppressWarnings(frm(o ~ x + cs(z), family = cumulative(),
                                data = d)))
out <- list()
for (nm in names(fits)) {
  for (pr in c("default", "flat")) {
    s <- suppressWarnings(suppressMessages(frm_sample(
      fits[[nm]], chains = 2, iter = 600, warmup = 300, seed = 1,
      cores = 1, refresh = 0,
      prior = if (pr == "flat") "flat")))
    out[[paste(nm, pr)]] <- s$draws
  }
}
saveRDS(out, sprintf("dev/optima-log/p1-nomo-%s.rds", arm))
cat("arm", arm, "saved", length(out), "draws matrices\n")
