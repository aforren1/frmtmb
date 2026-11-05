# Fix round 1, item 2, the rest: with a prior, SEs against quadrature and Laplace,
# profile bounds, update(), and the worst cost the refusal's new bound allows.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages(library(frmtmb))
sq <- function(e) suppressWarnings(suppressMessages(e))
ses <- function(f) sqrt(diag(frmtmb:::sdr_of(f)$cov.fixed))
set.seed(7)
G <- 30L
g <- factor(rep(seq_len(G), each = 5))
x <- rnorm(150)
y <- rbinom(150, 1, plogis(0.3 + 0.6 * x + rnorm(G, 0, 1.2)[g]))
d1 <- data.frame(y, x, g)
pri <- list(b = set_prior("normal(0, 0.1)", class = "b"),
            sd = set_prior("normal(0, 0.3)", class = "sd"),
            sd_gamma = set_prior("gamma(4, 2)", class = "sd"))
for (nm in names(pri)) {
  fi <- sq(frm(bf(y ~ x + (1 | g)), family = bernoulli(), data = d1, prior = pri[[nm]], importance = 1000L))
  fq <- sq(frm(bf(y ~ x + (1 | g)), family = bernoulli(), data = d1, prior = pri[[nm]], quadrature = TRUE))
  cat(sprintf("%-9s SE imp/quad %s  par z %s\n", nm, paste(format(ses(fi) / ses(fq), digits = 4), collapse = " "),
              paste(format((fi$opt$par - fq$opt$par) / ses(fq), digits = 3), collapse = " ")))
  ci <- sq(confint(fi, parm = "theta_1", method = "profile"))
  cq <- sq(confint(fq, parm = "theta_1", method = "profile"))
  cat(sprintf("          profile theta imp %.3f %.3f quad %.3f %.3f\n", ci[1, 1], ci[1, 2], cq[1, 1], cq[1, 2]))
  u <- sq(update(fi, importance = 200L))
  cat(sprintf("          update keeps the prior: objective at 200 draws %.3f, at 1000 %.3f\n", u$opt$objective, fi$opt$objective))
}
# the slope's lkj: gaussian SE against Laplace (Laplace is exact)
set.seed(5)
g2 <- factor(rep(seq_len(25L), each = 6))
x2 <- rnorm(150)
a <- rnorm(25, 0, 0.8); b <- rnorm(25, 0, 0.4)
y2 <- 0.3 + 0.5 * x2 + a[g2] + b[g2] * x2 + rnorm(150, 0, 0.5)
d2 <- data.frame(y = y2, x = x2, g = g2)
for (nm in c("lkj(4)", "lkj(0.5)")) {
  pr <- set_prior(nm, class = "cor")
  fi <- sq(frm(bf(y ~ x + (x | g)), family = gaussian(), data = d2, prior = pr, importance = 64L))
  fl <- sq(frm(bf(y ~ x + (x | g)), family = gaussian(), data = d2, prior = pr))
  cat(sprintf("gaussian %-9s obj-lap %+.3g maxpar %.3g SE imp/lap %s\n", nm, fi$opt$objective - fl$opt$objective,
              max(abs(fi$opt$par - fl$opt$par)), paste(format(ses(fi) / ses(fl), digits = 5), collapse = " ")))
}
