# Reviewer of lane illse, attack 6: item c's "at the maximum" judgment
# on fits stopped by an iteration limit (restarts = 0), against the true
# shortfall to a fully converged fit of the same model. Per fit: code,
# the warning's clause and stated gain, the true shortfall, whether the
# SE check then judged the fit (se_explained), and the SEs against the
# converged fit's.
#   REV_ARM=lane|base Rscript dev/illse-rev-conv.R <seeds>
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
a <- commandArgs(TRUE)
seeds <- if (length(a)) eval(parse(text = a[1])) else 501:503
models <- function(seed) {
  set.seed(seed)
  n <- 300
  x <- stats::rnorm(n); z <- stats::runif(n)
  g <- factor(rep(1:20, length.out = n))
  u <- stats::rnorm(20, 0, 0.6)[g]
  list(
    gauss = list(fo = y ~ x + z, fam = gaussian(),
                 d = data.frame(x, z, y = 1 + 2 * x - z + stats::rnorm(n))),
    pois_re = list(fo = y ~ x + (1 | g), fam = poisson(),
                   d = data.frame(x, g, y = stats::rpois(n, exp(0.5 + 0.4 * x + u)))),
    negbin = list(fo = y ~ x, fam = negbinomial(),
                  d = data.frame(x, y = stats::rnbinom(n, mu = exp(1 + 0.3 * x), size = 3))),
    student = list(fo = y ~ x, fam = student(),
                   d = data.frame(x, y = 1 + x + stats::rt(n, df = 4))),
    gamma_re = list(fo = y ~ x + (1 | g), fam = Gamma(link = "log"),
                    d = data.frame(x, g, y = stats::rgamma(n, 3, 3 / exp(0.2 + 0.3 * x + u)))))
}
num_after <- function(w, pat) {
  m1 <- regmatches(w, regexpr(paste0(pat, " [-0-9.e+]+"), w))
  if (!length(m1)) return(NA_real_)
  as.numeric(sub(".* ", "", m1))
}
for (seed in seeds) {
  ms <- models(seed)
  for (nm in names(ms)) {
    md <- ms[[nm]]
    full <- suppressWarnings(frm(md$fo, data = md$d, family = md$fam))
    f_best <- full$opt$objective
    se_full <- sqrt(diag(suppressWarnings(vcov(full, full = TRUE))))
    for (k in c(1, 2, 3, 4, 5, 6, 8, 10, 15)) {
      r <- cap(frm(md$fo, data = md$d, family = md$fam,
                   control = frmtmb_control(restarts = 0,
                                            optCtrl = list(iter.max = k,
                                                           eval.max = 1000))))
      if (inherits(r$v, "rev_err")) {
        cat(sprintf("%-8s s%d k %2d ERROR %s\n", nm, seed, k, substr(r$v, 1, 60)))
        next
      }
      f <- r$v
      w <- paste(r$w, collapse = " | ")
      clause <- if (grepl("so the estimates are at the maximum", w, fixed = TRUE)) "AT_MAX" else
        if (grepl("those of a point short of the maximum", w, fixed = TRUE)) "SHORT" else
          if (grepl("whether it is the maximum could not be measured", w, fixed = TRUE)) "UNUSABLE" else
            if (f$opt$convergence != 0) "none" else "code0"
      stated <- if (clause == "AT_MAX") num_after(w, "A Newton step from the stop gains") else
        if (clause == "SHORT") num_after(w, "would still gain") else NA
      true_short <- f$opt$objective - f_best
      se <- sqrt(diag(suppressWarnings(vcov(f, full = TRUE))))
      sr <- if (length(se) == length(se_full)) range(se / se_full) else c(NA, NA)
      other <- setdiff(r$w, grep("did not report convergence", r$w, value = TRUE))
      cat(sprintf(paste0("%-8s s%d k %2d | code %d | clause %-8s | stated %9.2e | ",
                         "true shortfall %9.2e | explained %-11s | lost %d | ",
                         "SE/full %.4f-%.4f | other warnings %d %s\n"),
                  nm, seed, k, f$opt$convergence, clause, stated, true_short,
                  format(f$cache$se_explained %||% "-"),
                  length(ns$sdr_of(f)$se_lost), sr[1], sr[2], length(other),
                  substr(paste(other, collapse = " / "), 1, 80)))
    }
  }
}
