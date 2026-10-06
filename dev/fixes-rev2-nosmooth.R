# Reviewer of lane fixes, re-check: models without a smooth, base
# against lane, bitwise (opt$par, objective, SEs). Writes one rds per lib.
#   Rscript dev/fixes-rev2-nosmooth.R <lib> <out.rds>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(5)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:30, 10)),
                w = runif(n, 1e3, 5e3))
d$y <- 1 + d$x + rnorm(30)[d$g] + rnorm(n)
d$c <- rpois(n, exp(0.3 + 0.4 * d$x + 0.3 * rnorm(30)[d$g]))
d$p <- 2 * exp(0.3 * d$x) + rnorm(n, 0, 0.2)
d$o <- cut(d$y, c(-Inf, 0, 1, 2, Inf), labels = FALSE)
d$yw <- 0.001 * d$w + d$x + rnorm(n)
ms <- list(
  lm = list(bf(y ~ x + z), gaussian()),
  lmm = list(bf(y ~ x + (1 + x | g)), gaussian()),
  glmm = list(bf(c ~ x + (1 | g)), poisson()),
  nb = list(bf(c ~ x + (1 | g)), negbinomial()),
  ls = list(bf(y ~ x, sigma ~ z), gaussian()),
  ord = list(bf(o ~ x + (1 | g)), cumulative()),
  nl = list(bf(p ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE), gaussian()),
  scaled = list(bf(yw ~ w + x), gaussian()))
out <- lapply(names(ms), function(k) {
  m <- ms[[k]]
  st <- if (k == "nl") list(beta = c(1, 0))
  f <- suppressWarnings(frm(m[[1]], data = d, family = m[[2]], start = st))
  list(par = f$opt$par, obj = f$opt$objective,
       se = suppressWarnings(fixef(f)[, "Est.Error"]), units = f$par_units)
})
names(out) <- names(ms)
saveRDS(out, a[2])
