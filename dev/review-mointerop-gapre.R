## Reviewer: (mo(x) | g) with an empty category, a converged fit.
##   Rscript dev/review-mointerop-gapre.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-mointerop/dev/review-mointerop-common.R")
pse <- function(fit, nd, ...) {
  as.data.frame(predictions(fit, newdata = nd, ...))[
    , c("estimate", "std.error")]
}
for (seed in c(51, 52)) {
  set.seed(seed)
  n <- 600
  G <- 20
  d <- data.frame(x = sample(c(1, 2, 5), n, TRUE), z = rnorm(n),
                  g = factor(sample(seq_len(G), n, TRUE)))
  u0 <- rnorm(G, 0, 0.6)
  u1 <- rnorm(G, 0, 0.3)
  d$y <- 1 + u0[d$g] + (0.4 + u1[d$g]) * (d$x - 1) + 0.5 * d$z + rnorm(n)
  fr <- show(paste("seed", seed, "(mo(x) | g), x in 1, 2, 5"),
             frm(y ~ mo(x) + z + (mo(x) | g), data = d))
  if (inherits(fr, "error")) next
  show("opt message", fr$opt$message)
  show("get_coef", get_coef(fr))
  show("get_vcov diag", diag(get_vcov(fr)))
  nd <- data.frame(x = c(1, 2, 3, 5), z = 0, g = factor(1, levels = 1:G))
  show("predictions conditional", pse(fr, nd))
  show("fitted conditional", fitted(fr, newdata = nd)[, 1:2])
  show("predictions conditional, observed rows only",
       pse(fr, nd[c(1, 2, 4), ]))
  show("predictions re_formula = NA, observed rows only",
       suppressWarnings(pse(fr, nd[c(1, 2, 4), ], re_formula = NA)))
  show("fitted re_formula = NA", fitted(fr, newdata = nd,
                                        re_formula = NA)[, 1:2])
  show("avg_slopes z", avg_slopes(fr, variables = "z"))
}
## The same without the random slope, to separate the two
set.seed(51)
n <- 600
G <- 20
d <- data.frame(x = sample(c(1, 2, 5), n, TRUE), z = rnorm(n),
                g = factor(sample(seq_len(G), n, TRUE)))
d$y <- 1 + rnorm(G, 0, 0.6)[d$g] + 0.4 * (d$x - 1) + 0.5 * d$z + rnorm(n)
fi <- show("(1 | g), x in 1, 2, 5", frm(y ~ mo(x) + z + (1 | g), data = d))
nd <- data.frame(x = c(1, 2, 3, 5), z = 0, g = factor(1, levels = 1:G))
show("predictions conditional", pse(fi, nd))
show("fitted conditional", fitted(fi, newdata = nd)[, 1:2])
show("predictions NA", suppressWarnings(pse(fi, nd, re_formula = NA)))
show("fitted NA", fitted(fi, newdata = nd, re_formula = NA)[, 1:2])
