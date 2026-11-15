# Review of lane wcheck2: sigma ~ 1 + fac on links other than log under
# a normal(0, 1) prior on the Intercept. The pair rule is decided on the
# log link only. With the identity link, sigma of cell b reaches 0 at a
# finite coefficient while the centered intercept and sigma of cell a
# stay finite; with softplus, a growing cell loses only log t.
#   Rscript dev/reviews/wcheck2-review-links.R <base|lane>
a <- commandArgs(TRUE)
lib_arg <- a[1]
source("dev/wcheck2-common.R")
set.seed(20261111)
d <- data.frame(obs = factor(seq_len(80)), x = rnorm(80),
                fac = factor(rep(c("a", "b"), 40)))
d$y <- 0.5 + 0.4 * d$x + rnorm(80)
fo <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + fac)
pic <- set_prior("normal(0, 1)", class = "Intercept", dpar = "sigma")
for (lk in c("log", "identity", "softplus")) {
  for (wb in c(1.2, 0.9)) {
    d$w <- ifelse(d$fac == "b", wb, 1)
    cat(sprintf("link %-8s w_b %.1f (E_b %+.0f): %s\n", lk, wb,
                40 * (wb - 1),
                outcome(frm(fo, data = d, prior = pic,
                            family = frmtmb:::fam_gaussian(link_sigma = lk)))))
  }
}
cat("DONE\n")
