# Reviewer check: brms 2.23.0's log density (rstan::log_prob on brms's
# own program, flat priors) at frmtmb's estimates, on cases the lane did
# not test:
#   hetgap  gaussian ar(cov = TRUE) + sigma ~ x2, unequal group lengths
#           AND interior time gaps (brms's normal_time_het_flex_lpdf)
#   stuhet  student arma(cov = TRUE) + sigma ~ x2, unequal lengths
#   stuar   student ar(cov = TRUE), constant sigma, gaps (hom_flex)
#   mare    gaussian ma(cov = TRUE) + (1 | h), joint density
#   mv      mvbf(y1 ~ x + ar(cov = TRUE), y2 ~ x), no rescor
# Usage: Rscript dev/review-sigmainn-brmslp.R <after|before>
arm <- commandArgs(TRUE)[1]
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sigmainn"
.libPaths(c(if (identical(arm, "after")) "C:/Users/adf44/source/r/wt-sigmainn-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = file.path(wt, "dev", "stan-cache"))
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
})
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
env <- new.env(parent = asNamespace("frmtmb"))
sys.source(file.path(wt, "tests", "testthat", "helper-brms.R"), envir = env)
# the helper reads a few names from the test file's scope
brms_lp_check <- env$brms_lp_check
mk <- function(seed, ng, nt, gaps = TRUE) {
  set.seed(seed)
  d <- do.call(rbind, lapply(seq_len(ng), function(i) {
    k <- sample(4:nt, 1)
    tt <- seq_len(nt)
    e <- as.numeric(arima.sim(list(ar = 0.6, ma = 0.3), nt, sd = 1))
    keep <- seq_len(k)
    if (gaps && i %% 3 == 0) keep <- setdiff(keep, 2L)
    data.frame(g = i, time = tt[keep], e = e[keep])
  }))
  d$g <- factor(d$g)
  d$x <- rnorm(nrow(d))
  d$x2 <- rnorm(nrow(d))
  d$h <- factor(sample(1:6, nrow(d), TRUE))
  d$y <- 1 + 0.7 * d$x + exp(-0.5 + 0.3 * d$x2) * d$e +
    c(-0.6, -0.3, 0, 0.2, 0.4, 0.5)[d$h]
  d$y1 <- d$y
  d$y2 <- 2 - 0.4 * d$x + rnorm(nrow(d))
  d
}
d <- mk(20261028, 30, 9)
cat("rows", nrow(d), "groups", nlevels(d$g), "lengths",
    paste(range(table(d$g)), collapse = "-"), "gapped groups",
    sum(tapply(d$time, d$g, function(t) any(diff(t) > 1))), "\n")
dn <- mk(20261029, 30, 9, gaps = FALSE)
cases <- list(
  hetgap = list(f = bf(y ~ x + ar(time, gr = g, cov = TRUE), sigma ~ x2),
                b = brms::bf(y ~ x + ar(time, gr = g, cov = TRUE), sigma ~ x2),
                fam = gaussian(), d = d, joint = FALSE),
  stuhet = list(f = bf(y ~ x + arma(time, gr = g, cov = TRUE), sigma ~ x2),
                b = brms::bf(y ~ x + arma(time, gr = g, cov = TRUE),
                             sigma ~ x2),
                fam = student(), bfam = brms::student(), d = dn, joint = FALSE),
  stuar = list(f = bf(y ~ x + ar(time, gr = g, cov = TRUE)),
               b = brms::bf(y ~ x + ar(time, gr = g, cov = TRUE)),
               fam = student(), bfam = brms::student(), d = dn, joint = FALSE),
  mare = list(f = bf(y ~ x + (1 | h) + ma(time, gr = g, cov = TRUE)),
              b = brms::bf(y ~ x + (1 | h) + ma(time, gr = g, cov = TRUE)),
              fam = gaussian(), d = dn, joint = TRUE),
  hetuneq = list(f = bf(y ~ x + ar(time, gr = g, cov = TRUE), sigma ~ x2),
                b = brms::bf(y ~ x + ar(time, gr = g, cov = TRUE), sigma ~ x2),
                fam = gaussian(), d = dn, joint = FALSE),
  mahet = list(f = bf(y ~ x + (1 | h) + ma(time, gr = g, cov = TRUE), sigma ~ x2),
               b = brms::bf(y ~ x + (1 | h) + ma(time, gr = g, cov = TRUE), sigma ~ x2),
               fam = gaussian(), d = dn, joint = TRUE))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  res <- tryCatch({
    fit <- frm(cs$f, family = cs$fam, data = cs$d)
    r <- NULL
    ok <- logical(0)
    withCallingHandlers(
      r <- brms_lp_check(cs$b, cs$bfam %||% cs$fam, cs$d, fit, joint = cs$joint),
      expectation = function(e) {
        ok <<- c(ok, inherits(e, "expectation_success"))
        invokeRestart("continue_test")
      })
    sprintf("lp %.10f ours %.10f diff %.3e max|grad| %.3e expectations %d/%d; sigma(fit) %s",
            r$lp, r$ours, r$lp - r$ours, r$max_grad, sum(ok), length(ok),
            paste(sprintf("%.6f", sigma(fit)), collapse = " "))
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-7s %s\n", nm, res))
}
