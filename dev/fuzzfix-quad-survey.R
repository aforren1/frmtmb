# Which multi-block quadrature fits leave random effects to Laplace, and
# which of them break summary()?  One row per arm.
#   Rscript dev/fuzzfix-quad-survey.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
dcum <- readRDS("dev/fuzzfix-out/data-cum.rds")

set.seed(20261020)
ng <- 24; nt <- 8; n <- ng * nt
d <- data.frame(g = factor(rep(seq_len(ng), each = nt)),
                tt = factor(rep(seq_len(nt), ng)),
                gb = factor(rep(c("i", "ii"), length.out = n)),
                gc = factor(rep(c("p", "p", "q", "q"), length.out = n)),
                x = rnorm(n), w = runif(n, 0.5, 2))
d$ga <- d$g
iab <- interaction(d$ga, d$gb, drop = TRUE)
iabc <- interaction(d$ga, d$gb, d$gc, drop = TRUE)
eta <- -0.2 + 0.4 * d$x + rnorm(ng, 0, 0.9)[d$g] +
  rnorm(nlevels(iab), 0, 0.9)[iab] + rnorm(nt, 0, 0.6)[d$tt]
d$yb <- rbinom(n, 1, plogis(eta))
d$yp <- rpois(n, exp(eta))
d$yo <- factor(cut(eta + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = FALSE),
               ordered = TRUE)
d$yg <- eta + rnorm(n, 0, exp(0.2 + rnorm(ng, 0, 0.4)[d$g]))

arms <- list(
  cum_fuzz = list(bf(y | weights(w) ~ 1 + x + (1 | ga/gb)) + cumulative(),
                  dcum),
  cum_fuzz_now = list(bf(y ~ 1 + x + (1 | ga/gb)) + cumulative(), dcum),
  cum_fuzz_ri = list(bf(y | weights(w) ~ 1 + x + (1 | ga)) + cumulative(),
                     dcum),
  bern_single = list(bf(yb ~ x + (1 | g)) + bernoulli(), d),
  bern_nested = list(bf(yb ~ x + (1 | ga/gb)) + bernoulli(), d),
  pois_nested = list(bf(yp ~ x + (1 | ga/gb)) + poisson(), d),
  pois_nested3 = list(bf(yp ~ x + (1 | ga/gb/gc)) + poisson(), d),
  pois_crossed = list(bf(yp ~ x + (1 | g) + (1 | tt)) + poisson(), d),
  cum_nested = list(bf(yo ~ x + (1 | ga/gb)) + cumulative(), d),
  sratio_nested = list(bf(yo ~ x + (1 | ga/gb)) + sratio(), d),
  acat_nested_w = list(bf(yo | weights(w) ~ x + (1 | ga/gb)) + acat(), d),
  grby3 = list(bf(yp ~ x + (1 | gr(g, by = gc))) + poisson(), d),
  gaus_mu_sigma = list(bf(yg ~ x + (1 | g), sigma ~ (1 | g)) + gaussian(),
                       d)
)
for (nm in names(arms)) {
  ar <- arms[[nm]]
  t0 <- proc.time()[[3]]
  fq <- tryCatch(suppressWarnings(frm(ar[[1]], data = ar[[2]],
                                      quadrature = TRUE)),
                 error = function(e) e)
  if (inherits(fq, "error")) {
    cat(sprintf("ARM %-14s fit-error: %s\n", nm,
                substr(conditionMessage(fq), 1, 150)))
    next
  }
  fl <- suppressWarnings(frm(ar[[1]], data = ar[[2]]))
  sm <- tryCatch(suppressWarnings(summary(fq)), error = function(e) e)
  vc <- tryCatch(suppressWarnings(vcov(fq)), error = function(e) e)
  cat(sprintf(paste0("ARM %-14s nb %3d left_random %3d logLik_quad %.8f ",
                     "logLik_lap %.8f diff %.3g summary %s vcov %s ",
                     "(%.1fs)\n"),
              nm, length(fq$estimates$b), length(fq$obj$env$random),
              as.numeric(logLik(fq)), as.numeric(logLik(fl)),
              as.numeric(logLik(fq)) - as.numeric(logLik(fl)),
              if (inherits(sm, "error")) paste0("ERROR[",
                conditionMessage(sm), "]") else "ok",
              if (inherits(vc, "error")) "ERROR" else "ok",
              proc.time()[[3]] - t0))
}
cat("END\n")
