# Review of lane nlstart: `b1 * exp(b2 * xb)`, xb in 0 to 3000, under
# cumulative() (review-nlstart-fresh.R body `bigx`). The probe's patterns
# sit at magnitude 1, where exp(b2 * xb) overflows. From starts of the
# four sign patterns at b2 = +-1e-3 (one unit of the body's scale) instead.
# Usage: REV_ARM=after Rscript dev/review-nlstart-bigx.R
source("dev/review-nlstart-common.R")
src <- readLines("dev/review-nlstart-fresh.R")
eval(parse(text = src[grep("^F <- list", src):(grep("^cat\\(\"ROWHEAD", src) - 1L)]))
for (seed in c(11, 22, 33, 44, 55)) {
  d <- gen("cumulative", "bigx", seed * 1000 + match("bigx", names(B)))
  fm <- form("cumulative", "bigx")
  f0 <- tt(frm(fm, family = cumulative(), data = d))
  res <- c()
  for (s1 in c(1, -1)) for (s2 in c(1, -1)) {
    f <- tt(frm(fm, family = cumulative(), data = d,
                start = list(beta = c("b1_(Intercept)" = s1,
                                      "b2_(Intercept)" = s2 * 1e-3))))
    res <- c(res, sprintf("(%+d,%+d) %.4f", s1, s2, llv(f)))
  }
  cat(sprintf("seed %d default %.4f conv %s nwarn %d multi %s | %s\n", seed,
              llv(f0), f0$opt$convergence, length(attr(f0, "w")),
              any(attr(f0, "mc") == "frmtmb_nl_multistart"),
              paste(res, collapse = "  ")))
}
