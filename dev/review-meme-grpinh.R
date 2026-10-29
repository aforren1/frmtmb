# Review of lane meme: a class-wide corme prior on a grouped correlated
# set, whose group row inherits it; are the prior draws that prior?
.libPaths(c("C:/Users/adf44/source/r/wt-meme-lib",
            "C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
source("C:/Users/adf44/source/r/frmtmb-wt-meme/dev/review-meme-shapes.R")
d <- rv_data()
f <- bf(y ~ me(xg, sg, gr = g) + me(wg, sg, gr = g))
run <- function(pr) {
  suppressWarnings(suppressMessages(frm_sample(
    f, data = d, family = gaussian(), prior = pr, sample_prior = "yes",
    chains = 2, iter = 1000, seed = 8, refresh = 0)))
}
ks <- function(v, a) suppressWarnings(stats::ks.test((v + 1) / 2, "pbeta",
                                                     a, a)$p.value)
ds <- run(set_prior("lkj(4)", class = "corme"))
pd <- prior_draws(ds)
cat("class lkj(4): names", names(pd), "\n")
if (!is.null(pd$Lme_1)) {
  cat(sprintf("  KS vs Beta(4,4): p = %.3g; vs Beta(1,1): p = %.3g\n",
              ks(pd$Lme_1, 4), ks(pd$Lme_1, 1)))
}
ts <- default_prior(f, data = d, route = "sample")
ts$prior[ts$class == "corme" & !nzchar(ts$group)] <- "lkj(4)"
print(as.data.frame(ts)[ts$class == "corme", c("prior", "class", "group")])
ds2 <- run(ts)
pd2 <- prior_draws(ds2)
cat("sample table, class row lkj(4), group row lkj(1) as printed: names",
    names(pd2), "\n")
if (!is.null(pd2$Lme_1)) {
  cat(sprintf("  KS vs Beta(4,4): p = %.3g; vs Beta(1,1): p = %.3g\n",
              ks(pd2$Lme_1, 4), ks(pd2$Lme_1, 1)))
}
rs <- suppressMessages(frmtmb.sample:::sample_resolve_priors(
  frm(f, data = d, dry_run = "objective"), ts))
cat("  sampled corme eta:", vapply(Filter(function(e) e$comp == "thetame",
                                          rs$ri$entries),
                                   function(e) e$dist$eta, 0), "\n")
cat("DONE\n")
