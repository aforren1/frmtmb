# Item 6 calibration: at which covariate scales does the DEFAULT fit fall
# short of autoscale = TRUE and of the reference (glm or the unscaled
# fit mapped back)? Run on the released build.
#   PREDFIX_ARM=base Rscript dev/predfix-scalescan.R > dev/predfix-log/scalescan-base.txt
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
ll <- function(f) if (inherits(f, "try-error")) NA_real_ else
  as.numeric(logLik(f))
scales <- 10^seq(-10, 10, by = 1)
mk <- function(seed, n = 250) {
  set.seed(seed)
  x <- rnorm(n)
  z <- runif(n)
  data.frame(
    x = x, z = z, f = factor(sample(letters[1:3], n, TRUE)),
    pois = rpois(n, exp(0.5 + 0.4 * x)),
    bin = rbinom(n, 1, plogis(0.3 + 0.8 * x)),
    gau = 1 + 0.5 * x + rnorm(n),
    gam = rgamma(n, 2, 2 / exp(0.2 + 0.3 * x)),
    nb = rnbinom(n, mu = exp(0.5 + 0.4 * x), size = 3)
  )
}
cases <- list(
  list(y = "pois", fam = poisson(), g = stats::poisson()),
  list(y = "bin", fam = bernoulli(), g = stats::binomial()),
  list(y = "gau", fam = gaussian(), g = stats::gaussian()),
  list(y = "gam", fam = Gamma(link = "log"), g = stats::Gamma("log")),
  list(y = "nb", fam = negbinomial(), g = NULL)
)
rhs <- c("0 + xs", "1 + xs", "xs + f")
rows <- list()
for (seed in 1:3) {
  d <- mk(seed)
  for (cs in cases) for (r in rhs) {
    ref <- NA_real_
    for (s in scales) {
      d$xs <- d$x * s
      fo <- stats::as.formula(paste(cs$y, "~", r))
      f0 <- try(quiet(frm(bf(fo), family = cs$fam, data = d)), silent = TRUE)
      f1 <- try(quiet(frm(bf(fo), family = cs$fam, data = d,
                          control = frmtmb_control(autoscale = TRUE))),
                silent = TRUE)
      gl <- if (!is.null(cs$g)) {
        g <- try(quiet(glm(fo, family = cs$g, data = d)), silent = TRUE)
        if (inherits(g, "try-error")) NA_real_ else as.numeric(logLik(g))
      } else NA_real_
      rows[[length(rows) + 1L]] <- data.frame(
        seed = seed, y = cs$y, rhs = r, scale = s, default = ll(f0),
        auto = ll(f1), glm = gl,
        conv0 = if (inherits(f0, "try-error")) NA else f0$opt$convergence)
    }
  }
}
res <- do.call(rbind, rows)
# the reference is the best logLik any arm reached at scale 1 for that
# design: ML is invariant to rescaling a column, so every scale shares it
ref <- aggregate(cbind(default, auto) ~ seed + y + rhs,
                 data = res[res$scale == 1, ], FUN = max)
res$ref <- pmax(ref$default, ref$auto)[match(
  paste(res$seed, res$y, res$rhs), paste(ref$seed, ref$y, ref$rhs))]
res$def_short <- res$ref - res$default
res$auto_short <- res$ref - res$auto
saveRDS(res, "C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-log/scalescan.rds")
options(width = 200)
cat("\nfits where the default falls more than 1e-6 short of the scale-1 optimum\n")
bad <- res[is.na(res$def_short) | res$def_short > 1e-6, ]
print(bad[, c("seed", "y", "rhs", "scale", "default", "auto", "glm",
              "def_short", "auto_short", "conv0")], digits = 10)
cat("\ncount short by log10 scale (default, autoscale)\n")
print(table(log10(res$scale), default_short = is.na(res$def_short) |
              res$def_short > 1e-6))
print(table(log10(res$scale), auto_short = is.na(res$auto_short) |
              res$auto_short > 1e-6))
