# Reviewer: default_prior() class/group rows against brms::get_prior() for
# the new structures. Data seed 20261008. Output: dev/ordinal-rev-log-priors.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261008)
n <- 200
d <- data.frame(x = rnorm(n), h = factor(sample(c("p", "q"), n, TRUE)))
d$y <- sample(1:4, n, TRUE); d$y2 <- sample(1:3, n, TRUE)
show <- function(lab, f, bfm, fam, bfam) {
  a <- as.data.frame(default_prior(f, data = d, family = fam))
  b <- as.data.frame(brms::get_prior(bfm, data = d, family = bfam))
  key <- function(z) sort(unique(paste(z$class, z$coef, z$group, z$resp, z$dpar,
                                       z$lb, sep = "|")))
  cat("\n==", lab, "\n frmtmb:", key(a), "\n brms:  ", key(b), "\n")
}
show("cumulative equidistant", y ~ x, brms::bf(y ~ x),
     cumulative(threshold = "equidistant"), brms::cumulative(threshold = "equidistant"))
show("sratio equidistant thres(gr = h)", y | thres(gr = h) ~ x,
     brms::bf(y | thres(gr = h) ~ x), sratio(threshold = "equidistant"),
     brms::sratio(threshold = "equidistant"))
show("acat sum_to_zero", y ~ x, brms::bf(y ~ x), acat(threshold = "sum_to_zero"),
     brms::acat(threshold = "sum_to_zero"))
show("cumulative disc ~ x", bf(y ~ x, disc ~ x), brms::bf(y ~ x, disc ~ x),
     cumulative(), brms::cumulative())
show("mv equidistant", bf(y ~ x) + bf(y2 ~ x),
     brms::bf(y ~ x) + brms::bf(y2 ~ x) + brms::set_rescor(FALSE),
     cumulative(threshold = "equidistant"), brms::cumulative(threshold = "equidistant"))
