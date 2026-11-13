# Short check 2 of lane wcheck: gaussian (1 | obs), E 0.5,
# normal(0, 1) on sigma, link_sigma identity. The exact (quadrature)
# log posterior along sigma (other parameters at the MAP of a fit
# with sigma held at 0.5) grows like 0.5 * log(1 / sigma); with the
# identity link the sampler's measure is d sigma, so the posterior mass
# of [s / 10, s] is about exp(f(s)) * s, which falls by 10^0.5 per
# decade: the posterior is proper there. The lane's sample route calls
# it improper.
#   Rscript dev/reviews/wcheck-short2-mass.R <lib>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(a[1], base, ul))
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(20261111)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n))
d$y <- 0.5 + 0.4 * d$x + rnorm(n)
d$w <- 1 + 0.5 / n
fam <- brmsfamily("gaussian", link_sigma = "identity")
pr <- set_prior("normal(0, 1)", class = "sigma")
ob <- q(frm(bf(y | weights(w) ~ x + (1 | obs)), family = fam, data = d,
            prior = pr, dry_run = "objective"))
ob <- if (!is.null(ob$obj)) ob$obj else ob
p <- ob$par
k <- which(names(p) == "betad")
free <- setdiff(seq_along(p), k)
s <- 10^-(1:8)
f <- vapply(s, function(sg) {
  p[k] <- sg
  fr <- function(r) {q2 <- p; q2[free] <- r
    v <- tryCatch(ob$fn(q2), error = function(e) NA_real_)
    if (!is.finite(v)) 1e10 else v}
  op <- nlminb(p[free], fr)
  if (op$objective < 1e9) p[free] <<- op$par
  -op$objective
}, 0)
cat("sigma", format(s), "\n-fn  ", format(f, digits = 7),
    "\nslope per unit log(1/sigma)",
    format(diff(f) / diff(-log(s)), digits = 3),
    "\nlog10 of exp(f) * sigma (mass per decade, profiled)",
    format((f + log(s)) / log(10), digits = 4), "\n")
cat("sample-route verdict:", tryCatch({
  fr <- q(frm(bf(y | weights(w) ~ x + (1 | obs)), family = fam, data = d,
              dry_run = "frame"))
  frmtmb:::check_weights_unbounded(fr$spec, fr, pr, exact = FALSE,
                                   route = "sample"); "pass"},
  error = function(e) substr(conditionMessage(e), 1, 60)), "\n")
cat("DONE\n")
