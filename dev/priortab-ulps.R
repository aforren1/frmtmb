# The rounding-sensitive assertions of test-prior-table-brms.R over 20
# ulp-perturbed responses, y * (1 + 2e-15 * k), k = 1..20.
# Rscript dev/priortab-ulps.R > dev/priortab-out/ulps.txt
.libPaths(c("C:/Users/adf44/source/r/wt-priortab-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(1)
n <- 120
d0 <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 12)))
d0$y <- 1 + 0.3 * d0$z + rnorm(10)[d0$g] * 0.3 + rnorm(n, sd = 0.3)
q <- function(e) suppressWarnings(suppressMessages(e))
ok <- matrix(NA, 20, 5, dimnames = list(NULL, c("box_identical", "at_0.2",
  "interior_rel", "roundtrip_rel", "roundtrip_null")))
for (k in 1:20) {
  d <- d0
  d$y <- d$y * (1 + 2e-15 * k)
  d$y3 <- d$y + 0.5 * d$z
  a <- q(frm(y3 ~ z, data = d, prior = set_prior("uniform(0, 0.2)",
                                                class = "b", coef = "z")))
  b <- q(frm(y3 ~ z, data = d, prior = set_prior("uniform(0, 0.2)",
             class = "b", coef = "z", lb = 0, ub = 0.2)))
  ml <- q(frm(y3 ~ z, data = d))
  ai <- q(frm(y3 ~ z, data = d, prior = set_prior("uniform(0, 2)",
                                                 class = "b", coef = "z")))
  zof <- function(f) unname(fixef(f)["z", "Estimate"])
  r0 <- q(frm(y ~ z + (1 | g), data = d))
  r1 <- q(frm(y ~ z + (1 | g), data = d,
              prior = default_prior(y ~ z + (1 | g), data = d)))
  ok[k, ] <- c(identical(a$opt$par, b$opt$par), identical(zof(a), 0.2),
               abs(zof(ai) - zof(ml)) / abs(zof(ml)),
               abs(r0$opt$objective - r1$opt$objective) /
                 abs(r0$opt$objective),
               is.null(r1$prior))
}
print(ok)
cat("max interior_rel", max(ok[, 3]), " max roundtrip_rel", max(ok[, 4]),
    "\n")
cat("all identical:", all(ok[, c(1, 2, 5)] == 1), "\n")
