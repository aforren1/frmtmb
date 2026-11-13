# Review of lane arspat: se() beside ar(cov = TRUE) on skew_normal, before
# (ARSPAT_LIB="") and after, and what brms 2.23.0 writes for it.
source("dev/reviews/arspat-review-common.R")
set.seed(7)
d <- expand.grid(t = 1:5, g = factor(1:30))
d$x <- rnorm(nrow(d)); d$s <- runif(nrow(d), 0.2, 0.9)
d$y <- 1 + 0.6 * d$x + rnorm(nrow(d))
chk("frmtmb skew_normal se ar cov=TRUE", {
  f <- frm(bf(y | se(s, sigma = TRUE) ~ x + ar(t, g, cov = TRUE)),
           family = skew_normal(), data = d)
  print(compat_row <- frm_compat("ar()", "se()"))
  logLik(f)
})
chk("frmtmb skew_normal se ar cov=FALSE", {
  f <- frm(bf(y | se(s, sigma = TRUE) ~ x + ar(t, g)),
           family = skew_normal(), data = d)
  logLik(f)
})
chk("frmtmb skew_normal se cosy", {
  f <- frm(bf(y | se(s, sigma = TRUE) ~ x + cosy(t, g)),
           family = skew_normal(), data = d)
  logLik(f)
})
if (nzchar(Sys.getenv("ARSPAT_BRMS"))) {
  chk("brms skew_normal se ar cov=TRUE", {
    sc <- brms::make_stancode(brms::bf(y | se(s, sigma = TRUE) ~ x +
                                         ar(t, g, cov = TRUE)),
                              family = brms::skew_normal(), data = d)
    cat(grep("se|err|ar", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
    NULL
  })
  chk("brms skew_normal se cosy", {
    sc <- brms::make_stancode(brms::bf(y | se(s, sigma = TRUE) ~ x +
                                         cosy(t, g)),
                              family = brms::skew_normal(), data = d)
    cat(grep("se|err|cosy", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
    NULL
  })
}
