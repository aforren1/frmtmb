# Lane seflags, punch round 1: degree-8 raw polynomials that keep every
# SE (seed 4 with (1 | g), seeds 2 and 10 without): frmtmb's SEs over
# lmer()'s (ML) or lm()'s with the ML residual variance.
#   Rscript dev/seflags-deg8cmp.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
cat("lib", find.package("frmtmb"), "\n")
for (cs in list(c(4, 1), c(2, 0), c(10, 0), c(7, 1))) {
  s <- cs[1]; re <- cs[2] == 1
  set.seed(s)
  d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
  for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
  rhs <- paste0("x", 1:8, collapse = " + ")
  fo <- as.formula(paste("y ~", rhs, if (re) "+ (1 | g)"))
  f <- suppressWarnings(frm(fo, data = d))
  se <- suppressWarnings(fixef(f)[, "Est.Error"])
  ref <- if (re) {
    sqrt(diag(as.matrix(vcov(suppressWarnings(lmer(fo, data = d,
                                                   REML = FALSE))))))
  } else {
    m <- lm(fo, data = d)
    sqrt(diag(vcov(m)) * (200 - 9) / 200)
  }
  cat("  se", format(se, digits = 4), "| ref", format(ref, digits = 4), "
")
  cat(sprintf("seed %2d re %-5s code %d | SE/ref %s\n", s, re,
              f$opt$convergence, paste(format(se / ref, digits = 3),
                                       collapse = " ")))
}
