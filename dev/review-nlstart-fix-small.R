# Review of lane nlstart, short check of the fix round: a covariate at scale
# 1e-4 inside a curved body, where the coefficient's maximum is near 1e4.
# Bodies `b1 * exp(b2 * xs)` and `b1 * (1 - exp(-b2 * xs))` (the latter's
# maximum has both signs negative), xs = x * 1e-4. Default fit against the
# best of 16 random starts (b2 at scale 1e4) and the truth, se = FALSE.
# Usage: REV_ARM=after|before Rscript dev/review-nlstart-fix-small.R [ref]
do_ref <- isTRUE(commandArgs(trailingOnly = TRUE)[1] == "ref")
source("dev/review-nlstart-common.R")
F <- list(gaussian = gaussian(), poisson = poisson(), bernoulli = bernoulli(),
          cumulative = cumulative())
B <- list(exp = list(f = "b1 * exp(b2 * xs)", tr = c(1.5, 6000)),
          sat = list(f = "b1 * (1 - exp(-b2 * xs))", tr = c(-1.5, -6000)))
cat("ROWHEAD body fam seed ll conv nwarn multi ref est\n")
for (bn in names(B)) for (fam in names(F)) for (seed in 1:5) {
  set.seed(seed * 31 + 7)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3))
  d$xs <- d$x * 1e-4
  lat <- 1.5 * exp(0.6 * d$x)
  d$y <- switch(fam, gaussian = lat + rnorm(n, 0, 0.5),
                poisson = rpois(n, lat), bernoulli = rbinom(n, 1, plogis(lat - 2)),
                cumulative = ordered(cut(lat + rlogis(n),
                                         c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4)))
  fm <- bf(as.formula(paste("y ~", B[[bn]]$f)), b1 + b2 ~ 1, nl = TRUE)
  r <- tt(frm(fm, family = F[[fam]], data = d))
  best <- NA
  if (do_ref) {
    set.seed(seed)
    sts <- c(lapply(1:16, function(i) c(rnorm(1, 0, 2), rnorm(1, 0, 6000))),
             list(B[[bn]]$tr))
    best <- max(vapply(sts, function(st) llv(tt(frm(fm, family = F[[fam]],
      data = d, se = FALSE, start = list(beta = c("b1_(Intercept)" = st[1],
                                                  "b2_(Intercept)" = st[2]))))),
      0), na.rm = TRUE)
  }
  cat(sprintf("ROW %s %s %d %.6f %s %d %d %.6f %s\n", bn, fam, seed, llv(r),
              if (iserr(r)) NA else r$opt$convergence, length(attr(r, "w")),
              as.integer(any(attr(r, "mc") == "frmtmb_nl_multistart")), best,
              if (iserr(r)) "ERR" else paste(format(r$estimates$beta, digits = 4),
                                             collapse = ",")))
  flush.console()
}
cat("DONE\n")
