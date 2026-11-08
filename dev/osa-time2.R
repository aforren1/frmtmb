# Lane osa, fix round: the cost of residuals(type = "osa") with (1 | g)
# against the 0.79.1 call on the same objective. The design of
# dev/review-osa-time.R and -time2.R: seed 601, groups of 20, x ~ N(0, 1),
# random intercept sd 0.7.
#
#   Rscript dev/osa-time2.R <family> <groups>
#
# family: poisson, negbinomial, bernoulli or binomial (10 trials); groups
# 25 (500 rows) or 100 (2000 rows). One call per arm, interleaved: "cdf"
# (the default route, its rows in blocks of about sqrt(n)), "old" (the
# 0.79.1 call), "one block" (the default route over every row at once,
# option frmtmb.osa_block = Inf; at 500 rows only, since at 2000 it is the
# half hour the review measured), and "cdf" again as the control, which
# must read about 1.
args <- commandArgs(TRUE)
fam <- args[1]
G <- as.integer(args[2])
.libPaths(c("C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib: frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
old <- function(fit) {
  RTMB::oneStepPredict(fit$obj, observation.name = ".frm_obs",
                       method = "oneStepGeneric", discrete = TRUE,
                       range = c(0, Inf), trace = FALSE)$residual
}
tm <- function(f) {
  t0 <- proc.time()[["elapsed"]]
  suppressWarnings(f())
  proc.time()[["elapsed"]] - t0
}
set.seed(601)
n <- G * 20L
d <- data.frame(g = factor(rep(seq_len(G), each = 20)), x = rnorm(n),
                k = 10L)
e <- rnorm(G, 0, 0.7)[d$g]
d$y <- switch(fam,
  poisson = rpois(n, exp(0.5 + 0.4 * d$x + e)),
  negbinomial = rnbinom(n, size = 3, mu = exp(0.5 + 0.4 * d$x + e)),
  bernoulli = rbinom(n, 1, plogis(0.5 + 0.4 * d$x + e)),
  binomial = rbinom(n, 10, plogis(0.5 + 0.4 * d$x + e)))
famobj <- frmtmb:::family_registry[[fam]]()
form <- if (fam == "binomial") {
  y | trials(k) ~ x + (1 | g)
} else y ~ x + (1 | g)
f1 <- suppressWarnings(frm(form, family = famobj, data = d))
t_new <- tm(function() residuals(f1, type = "osa"))
t_old <- tm(function() old(f1))
t_one <- if (G <= 25L) tm(function() {
  op <- options(frmtmb.osa_block = Inf)
  on.exit(options(op))
  residuals(f1, type = "osa")
}) else NA
t_ctl <- tm(function() residuals(f1, type = "osa"))
cat("BEGIN osa-time2\n")
cat(sprintf(paste("%-12s rows %4d (1 | g): cdf %7.1fs old %7.1fs",
                  "one block %7.1fs control %7.1fs | cdf/old %.2f",
                  "one block/old %.2f control/cdf %.2f\n"),
            fam, n, t_new, t_old, t_one, t_ctl, t_new / t_old,
            t_one / t_old, t_ctl / t_new))
cat("END osa-time2\n")
