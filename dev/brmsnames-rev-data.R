## Reviewer: shared data and model list for the independent brms checks.
## Data seed 2026. Sourced by brmsnames-rev-brms.R and -frm.R.
set.seed(2026)
n <- 360
rev_dd <- data.frame(
  x = stats::rnorm(n), z = stats::rnorm(n), w = stats::rnorm(n),
  sigma_z = stats::rnorm(n),
  f = factor(sample(c("a b", "c-d", "e"), n, TRUE)),
  g = factor(rep(1:12, 30)),
  h = factor(rep(c("lvl 1", "lvl-2", "lvl3"), 120)),
  h2 = factor(rep(c("p", "q", "r", "s", "t"), 72))
)
u <- matrix(stats::rnorm(12 * 4, 0, 0.5), 12)
rev_dd$y <- with(rev_dd, 1 + 0.5 * x + 0.2 * x^2 + u[g, 1] + u[g, 2] * x +
                   u[g, 3] * z + u[g, 4] * w + stats::rnorm(n, 0, 0.8))
rev_dd$y_a <- rev_dd$y
rev_dd$y2 <- with(rev_dd, 0.3 * x + u[g, 1] + stats::rnorm(n))
rev_dd$yn <- with(rev_dd, 2 * exp(-0.7 * abs(z)) * exp(u[g, 1] / 5) +
                    stats::rnorm(n, 0, 0.1))
rev_dd$az <- abs(rev_dd$z)

rev_models <- function(pkg) {
  bf <- get("bf", asNamespace(pkg))
  mvbf <- get("mvbf", asNamespace(pkg))
  list(
    C1 = list(f = bf(y ~ x + I(x^2) + f + x:f + (1 + x + z | g) + (1 | h))),
    C2 = list(f = bf(y ~ x + (1 + x + z + w | g))),
    C3 = list(f = bf(y ~ x + (1 | ID | g), sigma ~ z + (1 | ID | g))),
    C4 = list(f = mvbf(bf(y_a ~ x + (1 | p | g)), bf(y2 ~ x + (1 | p | g)),
                       rescor = FALSE)),
    C5 = list(f = bf(yn ~ a * exp(-b * az), a ~ 1 + (1 | g), b ~ 1,
                     nl = TRUE)),
    C6 = list(f = bf(y ~ sigma_z + (1 | g), sigma ~ z)),
    C7 = list(f = bf(y ~ x + (1 | g) + (0 + x | g) + (1 | g:h2))),
    C8 = list(f = bf(y ~ x + (1 | g)))
  )
}
