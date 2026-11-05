# R1: which fix stops the unbounded objective. Arms: none (no cap), NaN at a cap on
# the tape, or a box on the dispersion's betad at a cap. The review's designs
# (dev/review-latres-unb.R), seeds 1 to 30. The reference is the fit of the
# dispersion's limit on the same data, where the family has one: a beta_binomial
# without overdispersion is binomial. Usage:
#   Rscript dev/latres-r1.R <arm> <cap> <family> <outfile>
.libPaths(c("C:/Users/adf44/source/r/wt-latres-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
args <- commandArgs(TRUE)
arm <- args[[1]]; cap <- as.numeric(args[[2]]); famn <- args[[3]]; outf <- args[[4]]
options(warn = -1)
assignInNamespace("latres_disp_cap", if (arm == "nan") cap else Inf, "frmtmb")
ds <- list(
  bb = list(function(d) rbinom(nrow(d), 10, plogis(.3 * d$x + d$e)),
            y | trials(tr) ~ x + ar(t, g, cov = TRUE), beta_binomial(), binomial()),
  gamma = list(function(d) rgamma(nrow(d), 1000, 1000 / exp(.5 + .3 * d$x + d$e)),
               y ~ x + ar(t, g, cov = TRUE), Gamma(link = "log"), NULL),
  beta = list(function(d) {
    m <- plogis(.3 * d$x + d$e)
    pmin(pmax(rbeta(nrow(d), 1000 * m, 1000 * (1 - m)), 1e-4), 1 - 1e-4)
  }, y ~ x + ar(t, g, cov = TRUE), Beta(), NULL),
  nb = list(function(d) rpois(nrow(d), exp(1 + .3 * d$x + d$e)),
            y ~ x + ar(t, g, cov = TRUE), negbinomial(), poisson()),
  weibull = list(function(d) rweibull(nrow(d), 50, exp(.5 + .3 * d$x + d$e)),
                 y ~ x + ar(t, g, cov = TRUE), weibull(), NULL),
  zoib = list(function(d) {
    m <- plogis(.3 * d$x + d$e)
    v <- rbeta(nrow(d), 1000 * m, 1000 * (1 - m))
    u <- runif(nrow(d))
    ifelse(u < .1, rbinom(nrow(d), 1, .5), pmin(pmax(v, 1e-4), 1 - 1e-4))
  }, y ~ x + ar(t, g, cov = TRUE), zero_one_inflated_beta(), NULL))
D <- ds[[famn]]
first <- TRUE
for (s in 1:30) {
  set.seed(s)
  d <- expand.grid(t = 1:6, g = factor(1:25))
  d <- d[order(d$g, d$t), ]
  n <- nrow(d)
  d$x <- rnorm(n)
  d$e <- as.vector(replicate(25, arima.sim(list(ar = .5), 6, sd = .5)))
  d$tr <- 10L
  d$y <- D[[1]](d)
  w <- character(0)
  pr <- if (arm == "box") set_prior("", class = if (famn %in% c("bb", "beta", "zoib")) "phi" else "shape", ub = exp(cap))
  f <- withCallingHandlers(
    tryCatch(frm(D[[2]], data = d, family = D[[3]], prior = pr),
             error = function(e) e),
    warning = function(c) {
      w <<- c(w, conditionMessage(c))
      invokeRestart("muffleWarning")
    }, message = function(c) invokeRestart("muffleMessage"))
  out <- data.frame(arm = arm, cap = cap, fam = famn, seed = s, err = "", ll = NA_real_,
                    code = NA, nw = length(w), ldisp = NA_real_, sderr = NA_real_,
                    ar = NA_real_, ref = NA_real_, ref_s = NA_real_)
  if (inherits(f, "error")) {
    out$err <- substr(conditionMessage(f), 1, 60)
  } else {
    par <- f$opt$par
    out$ll <- as.numeric(logLik(f))
    out$code <- f$opt$convergence
    out$ldisp <- if (any(names(par) == "betad")) par[names(par) == "betad"][1] else NA
    th <- f$estimates$theta
    out$sderr <- exp(th[1]); out$ar <- th[2] / sqrt(1 + th[2]^2)
    if (!is.null(D[[4]])) {
      r <- tryCatch(suppressMessages(frm(D[[2]], data = d, family = D[[4]])),
                    error = function(e) NULL)
      if (!is.null(r)) {
        out$ref <- as.numeric(logLik(r)); out$ref_s <- exp(r$estimates$theta[1])
      }
    }
  }
  write.table(out, outf, sep = "\t", row.names = FALSE, quote = FALSE,
              col.names = first, append = !first)
  first <- FALSE
}
