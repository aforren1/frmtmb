# Reviewer, fix round 1: which shapes does quad_check_blocks() refuse,
# and with which message? One arm per process.
#   Rscript dev/fuzzfix-rev-r1-quad.R <lib|base> <arm>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
arm <- a[2]
set.seed(31)
ng <- 12; m <- 4
d <- expand.grid(r = seq_len(m), gb = factor(c("i", "ii")),
                 ga = factor(seq_len(ng)))
n <- nrow(d)
d$x <- rnorm(n)
d$xs <- runif(n, 0, 10)
d$tt <- factor(rep(1:3, length.out = n))
d$fA <- factor(ifelse(as.integer(d$ga) %% 2 == 0, "even", "odd"))
cell <- interaction(d$ga, d$gb, drop = TRUE)
eta <- 0.3 + 0.5 * d$x + rnorm(ng, 0, 0.8)[d$ga] +
  rnorm(nlevels(cell), 0, 0.6)[cell]
d$y <- rbinom(n, 1, plogis(eta))
d$k <- rpois(n, exp(eta))
d$yg <- eta + rnorm(n, 0, 0.7)
B <- bernoulli()
arms <- list(
  single = list(bf(y ~ x + (1 | ga)), B),
  single_pois = list(bf(k ~ x + (1 | ga)), poisson()),
  wide = list(bf(y ~ x + (1 + x | ga)), B),
  dbar = list(bf(y ~ x + (x || ga)), B),
  smooth_ri = list(bf(y ~ s(xs) + (1 | ga)), B),
  by_split = list(bf(y ~ x + (1 | gr(ga, by = fA))), B),
  mv = list(mvbf(bf(y ~ x + (1 | ga)) + B, bf(k ~ x + (1 | ga)) + poisson()),
            NULL),
  mv_two_groups = list(mvbf(bf(y ~ x + (1 | ga)) + B,
                            bf(k ~ x + (1 | tt)) + poisson()), NULL),
  nested = list(bf(y ~ x + (1 | ga/gb)), B),
  crossed = list(bf(y ~ x + (1 | ga) + (1 | tt)), B),
  mu_sigma = list(bf(yg ~ x + (1 | ga), sigma ~ (1 | ga)), gaussian()),
  mu_sigma_other = list(bf(yg ~ x + (1 | ga), sigma ~ (1 | tt)), gaussian())
)
ar <- arms[[arm]]
t0 <- proc.time()[[3]]
f <- tryCatch(suppressWarnings(suppressMessages(
  frm(ar[[1]], family = ar[[2]], data = d, quadrature = TRUE))),
  error = function(e) e)
el <- proc.time()[[3]] - t0
if (inherits(f, "error")) {
  cat(sprintf("ARM %s REFUSED (%.1fs): %s\n", arm, el, conditionMessage(f)))
} else {
  sm <- tryCatch(suppressWarnings(suppressMessages(summary(f))),
                 error = function(e) e)
  cat(sprintf("ARM %s ok (%.1fs) blocks %d left_random %d logLik %.10f summary %s\n",
              arm, el, length(f$frame$re_blocks), length(f$obj$env$random),
              as.numeric(logLik(f)),
              if (inherits(sm, "error")) conditionMessage(sm) else "ok"))
}
