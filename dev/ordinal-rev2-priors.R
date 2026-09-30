# Reviewer: default_prior() class/group rows against brms::get_prior() for
# the new structures. Data seed 20261008. Output: dev/ordinal-rev2-log-priors.txt (re-check copy)
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
show("cumulative equidistant thres(gr = h)", y | thres(gr = h) ~ x,
     brms::bf(y | thres(gr = h) ~ x), cumulative(threshold = "equidistant"),
     brms::cumulative(threshold = "equidistant"))
show("hurdle equidistant", y ~ x, brms::bf(y ~ x),
     hurdle_cumulative(threshold = "equidistant"),
     brms::hurdle_cumulative(threshold = "equidistant"))
show("cumulative flexible thres(gr = h)", y | thres(gr = h) ~ x,
     brms::bf(y | thres(gr = h) ~ x), cumulative(), brms::cumulative())
show("acat sum_to_zero thres(gr = h)", y | thres(gr = h) ~ x,
     brms::bf(y | thres(gr = h) ~ x), acat(threshold = "sum_to_zero"),
     brms::acat(threshold = "sum_to_zero"))
# a per-level Intercept row is usable: the prior moves only that level
f0 <- frm(y | thres(gr = h) ~ x, family = cumulative(), data = d)
f1 <- frm(y | thres(gr = h) ~ x, family = cumulative(), data = d,
          prior = set_prior("normal(5, 0.01)", class = "Intercept", group = "p"))
th <- function(f) frmtmb:::ord_threshold_values(family(f), f$estimates$tau_raw)
cat("\nper-level Intercept prior (group p at N(5, .01)): thresholds\n")
print(rbind(ML = th(f0), prior_p = th(f1)))
