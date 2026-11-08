# Review of lane nlstart, short check: (1) the rows of fix-order5.log where
# the order of the terms changed the fit: each spelling with the probe on
# and off, against the best of 16 random starts; (2) the sat body at
# covariate scale 1e-4 (fix-small.R): the probe's start sizes, and the four
# sign patterns from b1 = +-1, b2 = +-2e4 as starts.
# Usage: REV_ARM=after Rscript dev/review-nlstart-fix-detail.R
source("dev/review-nlstart-common.R")
src <- readLines("dev/review-nlstart-fix-order5.R")
eval(parse(text = src[grep("^F <- list", src):(grep("^cat\\(\"ROWHEAD", src) - 1L)]))
rows <- list(c("six", "cumulative", 503), c("tie", "bernoulli", 502),
             c("tie", "cumulative", 502), c("tie", "cumulative", 503),
             c("tie", "poisson", 503))
for (r in rows) {
  bn <- r[1]; fam <- r[2]; seed <- as.numeric(r[3])
  terms <- switch(bn, six = terms6, tie = termsT)
  d <- gen(fam, seed, tie = bn == "tie")
  out <- c()
  best <- -Inf
  for (rv in c(FALSE, TRUE)) {
    s <- spell(terms, 6, maps$b, rv)
    fm <- bf(as.formula(paste("y ~", s$body)), as.formula(paste(s$decl, "~ 1")),
             nl = TRUE)
    on <- tt(frm(fm, family = F[[fam]], data = d))
    off <- withr::with_options(list(frmtmb.nl_sign_probe = FALSE),
                               tt(frm(fm, family = F[[fam]], data = d)))
    out <- c(out, sprintf("%s: on %.4f (w%d m%d) off %.4f", if (rv) "reversed" else "in order",
                          llv(on), length(attr(on, "w")),
                          as.integer(any(attr(on, "mc") == "frmtmb_nl_multistart")),
                          llv(off)))
    if (!rv) {
      set.seed(seed)
      for (i in 1:16) {
        st <- setNames(rnorm(6) * if (i > 8) 2 else 1, paste0("b", 1:6, "_(Intercept)"))
        best <- max(best, llv(tt(frm(fm, family = F[[fam]], data = d, se = FALSE,
                                     start = list(beta = st)))), na.rm = TRUE)
      }
    }
  }
  cat(sprintf("%s %s %d | %s | best of 16 starts %.4f\n", bn, fam, seed,
              paste(out, collapse = " | "), best))
}
cat("\n== sat body, xs = x * 1e-4, gaussian (fix-small.R data)\n")
for (seed in 1:3) {
  set.seed(seed * 31 + 7)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3))
  d$xs <- d$x * 1e-4
  d$y <- 1.5 * exp(0.6 * d$x) + rnorm(n, 0, 0.5)
  fm <- bf(y ~ b1 * (1 - exp(-b2 * xs)), b1 + b2 ~ 1, nl = TRUE)
  f0 <- tt(frm(fm, data = d))
  fr <- f0$frame
  lp <- Filter(function(l) !is.null(l[["nl_body"]]), fr[["linpreds"]])[[1]]
  probed <- list(list(par = "beta", at = 1L, j = 1L),
                 list(par = "beta", at = 2L, j = 2L))
  sc <- frmtmb:::nl_probe_scale(fr, lp, f0$estimates, probed)
  res <- c()
  for (s1 in c(1, -1)) for (s2 in c(1, -1)) {
    f <- tt(frm(fm, data = d, start = list(beta = c("b1_(Intercept)" = s1,
                                                    "b2_(Intercept)" = s2 * 2e4))))
    res <- c(res, sprintf("(%+d,%+d) %.4f", s1, s2, llv(f)))
  }
  cat(sprintf("seed %d default %.4f conv %s nwarn %d | probe sizes b1 %s b2 %s | %s\n",
              seed, llv(f0), f0$opt$convergence, length(attr(f0, "w")),
              format(sc["mag", 1]), format(sc["mag", 2]), paste(res, collapse = "  ")))
}
