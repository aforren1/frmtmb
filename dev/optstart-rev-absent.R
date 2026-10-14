# Reviewer of lane optstart: is the start 0.70.0's, bit for bit, where
# the rule does not apply but the unit changed (weights, mi(),
# zero_inflated_asym_laplace)? Prints the cold start and the predictor
# units in hex; run on both libraries and diff.
#   Rscript dev/optstart-rev-absent.R <lib> <label>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(3)
n <- 300
d <- data.frame(x = rnorm(n), g = factor(rep(1:15, 20)),
                w = rep(c(0.2, 1, 5), 100), cn = rbinom(n, 1, 0.2))
d$y <- 20 + 30 * d$x + 50 * (d$w > 1) + rnorm(15, 0, 3)[d$g] +
  rnorm(n, 0, 1 / sqrt(d$w))
d$ym <- d$y
d$ym[sample(n, 80)] <- NA
d$yz <- ifelse(rbinom(n, 1, 0.3) == 1, 0, d$y)
cases <- list(
  "weights + s(x), rule off" = list(bf(y | weights(w) ~ s(x) + (1 | g)),
                                    gaussian()),
  "weights + cens, rule off" = list(bf(y | weights(w) + cens(cn) ~ x +
                                         (1 | g)), gaussian()),
  "mi() + s(x), rule off" = list(bf(ym | mi() ~ s(x) + (1 | g)), gaussian()),
  "mi() intercept only, rule off" = list(bf(ym | mi() ~ 1 + (1 | g)),
                                         gaussian()),
  "zi_asym_laplace intercept only, rule off" = list(
    bf(yz ~ 1 + (1 | g), quantile = 0.5), zero_inflated_asym_laplace()))
hex <- function(v) paste(sprintf("%a", unname(v)), collapse = " ")
for (nm in names(cases)) {
  cs <- cases[[nm]]
  st <- tryCatch(unclass(par_template(cs[[1]], data = d, family = cs[[2]])),
                 error = function(e) NULL)
  fr <- frm(cs[[1]], data = d, family = cs[[2]], dry_run = "frame")
  u <- ns$predictor_units(fr)
  cat(sprintf("%s\n  betad %s\n  theta %s\n  units %s\n", nm,
              hex(st$betad), hex(st$theta), hex(u)))
}
