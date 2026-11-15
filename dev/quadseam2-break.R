# Item 561: does single-block quadrature still break down on a variance
# component near zero? Fits with a true group sd of 0 (and 0.05) over
# families whose inverse link exponentiates, quadrature = TRUE, counting
# the refusal "could not marginalize this model" and other errors.
#   Rscript dev/quadseam2-break.R <libpath|base> <seeds>
a <- commandArgs(TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 1:20
gen <- function(fam, s, sdg) {
  set.seed(s)
  G <- 15; m <- 4
  d <- data.frame(g = factor(rep(seq_len(G), each = m)), x = rnorm(G * m))
  u <- rnorm(G, 0, sdg)[d$g]
  eta <- 0.3 + 0.5 * d$x + u
  d$y <- switch(fam,
    poisson = rpois(nrow(d), exp(eta)),
    Gamma = rgamma(nrow(d), 3, rate = 3 / exp(eta)),
    beta = { mu <- plogis(eta); rbeta(nrow(d), mu * 8, (1 - mu) * 8) },
    negbinomial = rnbinom(nrow(d), mu = exp(eta), size = 3),
    bernoulli = rbinom(nrow(d), 1, plogis(eta)),
    lognormal = rlnorm(nrow(d), eta, 0.5))
  d
}
fams <- list(poisson = poisson(), Gamma = Gamma(link = "log"),
             beta = Beta(), negbinomial = negbinomial(),
             bernoulli = bernoulli(), lognormal = lognormal())
`%||%` <- function(a, b) if (is.null(a)) b else a
tab <- list()
for (fn in names(fams)) for (sdg in c(0, 0.05)) for (s in seeds) {
  d <- gen(fn, s, sdg)
  r <- tryCatch({
    f <- suppressMessages(suppressWarnings(frm(y ~ x + (1 | g),
           family = fams[[fn]], data = d, quadrature = TRUE)))
    sprintf("ok sd %.3g", exp(f$estimates$theta[1]))
  }, error = function(e) {
    m <- conditionMessage(e)
    if (grepl("could not marginalize", m)) "BREAKDOWN" else
      paste("ERROR", substr(gsub("\\s+", " ", m), 1, 100))
  })
  cat(sprintf("ROW\t%s\t%.2f\t%d\t%s\n", fn, sdg, s, r))
  key <- paste(fn, sdg, sub(" .*", "", r))
  tab[[key]] <- (tab[[key]] %||% 0) + 1
}

for (k in names(tab)) cat("SUM", k, tab[[k]], "\n")
cat("END\n")
