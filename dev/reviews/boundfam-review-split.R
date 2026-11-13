# Review of lane boundfam: message/warning split on slopes of family
# parameters (sigma ~ g, hu ~ x) beside the zi slope of M1.
#   Rscript dev/reviews/boundfam-review-split.R <lane|base>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-boundfam-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-boundfam/dev/reviews/boundfam-review-fp.R")
eval(parse(text = src[grep("^conds <- ", src):(grep("^ok <- ", src) - 1L)]))
F <- function(...) conds(frm(...))
# hurdle probability separating in x: zeros only below 0
for (s in 1:2) {
  set.seed(s)
  d <- data.frame(x = rnorm(200))
  d$y <- ifelse(d$x < 0, 0L, 1L + rpois(200, 2))
  show(paste0("hu sep s", s), F(bf(y ~ 1, hu ~ x), family = hurdle_poisson(), data = d))
}
# zero inflation run off in one group of a factor (no excess zeros in b)
for (s in 1:2) {
  set.seed(s)
  d <- data.frame(g = factor(rep(c("a", "b"), each = 150)))
  z <- rbinom(300, 1, 0.3) * (d$g == "a")
  d$y <- ifelse(z == 1, 0L, rpois(300, 4))
  show(paste0("zi ~ g s", s), F(bf(y ~ 1, zi ~ g), family = zero_inflated_poisson(), data = d))
}
# negbinomial shape ~ g, one group poisson-like
for (s in 1:2) {
  set.seed(s)
  d <- data.frame(g = factor(rep(c("a", "b"), each = 200)))
  d$y <- ifelse(d$g == "a", MASS::rnegbin(400, 4, 2), rpois(400, 4))
  show(paste0("shape ~ g s", s), F(bf(y ~ 1, shape ~ g), family = negbinomial(), data = d))
}
