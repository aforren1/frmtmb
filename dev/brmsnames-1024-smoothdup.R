## A smooth's unpenalized column written in both spellings, s(x).fx1
## and brms's sx_1, with two different tight priors (review short check,
## R1). On 0.74.0 both are accepted and one is dropped without a word;
## on the lane the pair is refused as a duplicate. Data: the tests'
## bns_prior_data(), seed 31.
##   Rscript dev/brmsnames-1024-smoothdup.R base|lane > dev/brmsnames-1024-log/smoothdup-<arm>.txt
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-brmsnames-lib", base)
} else base)
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
cat("arm", arm, dirname(system.file(package = "frmtmb")), "\n")
set.seed(31)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), xp = runif(n, 0.5, 2),
                k = factor(sample(1:3, n, replace = TRUE)),
                g = factor(rep(1:10, length.out = n)))
d$y <- 1 + d$x + 0.3 * d$z^2 + 0.5 * log(d$xp) + rnorm(10)[d$g] + rnorm(n)
fs <- bf(y ~ s(x) + z)
p <- set_prior("normal(3, 0.001)", class = "b", coef = "s(x).fx1") +
  set_prior("normal(-3, 0.001)", class = "b", coef = "sx_1")
v <- tryCatch(q(validate_prior(p, fs, data = d)), error = function(e) e)
if (inherits(v, "error")) {
  cat("validate_prior:", conditionMessage(v), "\n")
} else {
  v <- as.data.frame(v)
  u <- v[v$source == "user", c("prior", "class", "coef")]
  cat("validate_prior accepts both; user rows:\n")
  print(u, row.names = FALSE)
}
f <- tryCatch(q(frm(fs, data = d, prior = p)), error = function(e) e)
if (inherits(f, "error")) {
  cat("frm:", conditionMessage(f), "\n")
} else {
  cat("frm accepts both; fixef sx_1 =",
      format(fixef(f)["sx_1", "Estimate"], digits = 6), "\n")
}
