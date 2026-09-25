# The caller's RNG stream around simulate() and pp_check() at newdata.
#   Rscript dev/simnewdata-review/rv-rng.R > .../log/rng.txt
source("dev/simnewdata-review/rv-prelude.R")
set.seed(7)
n <- 120
d <- data.frame(x = rnorm(n), g = factor(rep(1:12, 10)))
d$y <- 0.5 + 0.8 * d$x + rnorm(12)[d$g] + rnorm(n, 0, 0.5)
fit <- frm(bf(y ~ x + (1 + x | g)), data = d)
nd <- data.frame(x = c(0, 1, 2), g = factor(c("1", "new", "new")))
cases <- list(
  "seed, newdata, allow_new_levels" = function() simulate(
    fit, nsim = 5, seed = 1, newdata = nd, allow_new_levels = TRUE),
  "seed, NA, newdata" = function() simulate(fit, nsim = 5, seed = 1,
                                            re_formula = NA, newdata = nd),
  "seed, partial ~(1|g)" = function() simulate(fit, nsim = 5, seed = 1,
                                               re_formula = ~ (1 | g)),
  "seed, partial at newdata" = function() simulate(
    fit, nsim = 5, seed = 1, re_formula = ~ (1 | g), newdata = nd,
    allow_new_levels = TRUE),
  "seed, error mid-call" = function() try(simulate(
    fit, nsim = 5, seed = 1, newdata = data.frame(x = 1, g = "zz")),
    silent = TRUE)
)
for (nm in names(cases)) {
  set.seed(99)
  before <- .Random.seed
  cases[[nm]]()
  cat(sprintf("%-34s caller stream unchanged: %s\n", nm,
              identical(before, .Random.seed)))
}
# the same seed gives the same draws whatever the caller's state
set.seed(1); a <- simulate(fit, nsim = 3, seed = 5, newdata = nd,
                           allow_new_levels = TRUE)
set.seed(2); b <- simulate(fit, nsim = 3, seed = 5, newdata = nd,
                           allow_new_levels = TRUE)
cat("seeded draws independent of caller state:", identical(a, b), "\n")
# a row's draws do not depend on rows after it? (not promised; record)
a1 <- simulate(fit, nsim = 3, seed = 5, newdata = nd[1, ])
cat("row 1 alone equals row 1 of the 3-row call:",
    identical(unlist(a1[1, ]), unlist(a[1, ])), "\n")
# unseeded: how much of the stream does one replicate at newdata take,
# against one at the fitted rows (rnorm count via a counter)
count_draws <- function(expr) {
  set.seed(3); expr; s1 <- .Random.seed
  set.seed(3); k <- 0L
  repeat { if (identical(.Random.seed, s1)) break; runif(1); k <- k + 1L
    if (k > 1e6) return(NA) }
  k
}
cat("stream draws, simulate(nsim = 1) fitted rows:",
    count_draws(simulate(fit, nsim = 1)), "\n")
cat("stream draws, simulate(nsim = 1, newdata = nd, allow):",
    count_draws(simulate(fit, nsim = 1, newdata = nd,
                         allow_new_levels = TRUE)), "\n")
cat("stream draws, simulate(nsim = 1, re_formula = NA, newdata = nd[1,]):",
    count_draws(simulate(fit, nsim = 1, re_formula = NA,
                         newdata = nd[1, ])), "\n")
cat("stream draws, predict(newdata = nd[1, ]) conditional:",
    count_draws(predict(fit, newdata = nd[1, ], ndraws = 1)), "\n")
