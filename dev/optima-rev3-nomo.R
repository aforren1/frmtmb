# Reviewer of lane optima, final check: check_laplace() on a fit with no
# mo() term, base against lane. Writes the returned table and the
# messages to an rds per arm and run; dev/optima-rev3-nomo-cmp.R
# compares them. Each run is pinned to one core (pitfall 21).
#   Rscript dev/optima-rev3-nomo.R base|lane <tag>
arm <- commandArgs(TRUE)[1]
tag <- commandArgs(TRUE)[2]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
out <- list()
set.seed(4)
dd <- data.frame(x = rnorm(120), g = factor(rep(1:30, 4)))
dd$y <- rbinom(120, 1, plogis(0.3 + 0.5 * dd$x + rnorm(30, 0, 1)[dd$g]))
fits <- list(
  glmm = frm(bf(y ~ x + (1 | g)) + frmtmb::bernoulli(), data = dd),
  map = frm(bf(y ~ x + (1 | g)) + frmtmb::bernoulli(), data = dd,
            prior = set_prior("normal(0, 1)", class = "b")))
for (nm in names(fits)) {
  msg <- character()
  r <- withCallingHandlers(
    check_laplace(fits[[nm]], chains = 2, iter = 600, seed = 3,
                  refresh = 0, cores = 1),
    message = function(m) {
      msg <<- c(msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) invokeRestart("muffleWarning"))
  out[[nm]] <- list(tab = r, msg = msg)
}
saveRDS(out, sprintf(
  "C:/Users/adf44/source/r/frmtmb-wt-optima/dev/optima-rev3-out/nomo-%s-%s.rds",
  arm, tag))
cat("done", arm, tag, find.package("frmtmb.sample"), "\n")
