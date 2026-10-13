# Reviewer of lane seflags, claim 4's joining rule on identified
# directions: a model with one exact alias (a truly flat direction) and
# an unrelated near-collinear but identified pair z2 = z1 + N(0, s). Does
# the pair lose its SEs as "flat" (the data do not determine them)? The
# reference is the same likelihood without the alias.
#   Rscript dev/seflags-rev2-join.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev2-common.R")
mk <- function(s, seed = 1) {
  set.seed(seed)
  n <- 200
  x <- stats::rnorm(n) + 0.1
  z1 <- stats::rnorm(n)
  d <- data.frame(x = x, z1 = z1, z2 = z1 + stats::rnorm(n, 0, s),
                  g = factor(rep(1:10, n / 10)))
  d$y <- 1 + 0.5 * d$x + sin(2 * d$x) + 0.3 * d$z1 + 0.3 * d$z2 +
    stats::rnorm(10, 0, 0.5)[d$g] + stats::rnorm(n, 0, 0.3)
  d
}
row <- function(lab, r) {
  fit <- r$v
  se <- se_outer(fit)
  zz <- se[names(se) %in% c("z1", "z2", "a_z1", "a_z2")]
  cat(sprintf("%-34s code %d | SE z1 z2 %-22s | warn %d | lost %s\n", lab,
              fit$opt$convergence, paste(format(zz, digits = 6),
                                         collapse = " "),
              length(r$w), lost_of(fit)))
}
for (seed in 1:3) for (s in c(1e-4, 3e-5, 1e-5, 3e-6, 1e-6)) {
  d <- mk(s, seed)
  tag <- sprintf("seed %d s %g", seed, s)
  row(paste(tag, "alias"), cap(frm(y ~ x + s(x) + z1 + z2 + (1 | g),
                                   data = d)))
  row(paste(tag, "no alias"), cap(frm(y ~ s(x) + z1 + z2 + (1 | g),
                                      data = d)))
  row(paste(tag, "nl a+b"), cap(frm(bf(y ~ a + b, a ~ 1 + x + z1 + z2 +
                                         (1 | g), b ~ 1, nl = TRUE),
                                    data = d)))
  m <- suppressWarnings(lme4::lmer(y ~ x + z1 + z2 + (1 | g), data = d,
                                   REML = FALSE))
  cat(sprintf("%-34s lmer y ~ x + z1 + z2 + (1|g): SE z1 z2 %s\n", tag,
              paste(format(sqrt(diag(as.matrix(vcov(m))))[3:4], digits = 6),
                    collapse = " ")))
}
