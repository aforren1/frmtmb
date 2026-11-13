# One-block quadrature accuracy against glmer(nAGQ = 25), base and lane:
# bernoulli (groups of 4, sd 1.2, the design of test-quadrature-defects.R
# "modes match glmer"), binomial of 10 trials and poisson (60 groups of
# 6, sd 0.3), seeds 1..S. Per fit: logLik - glmer25, the re-anchorings
# made (lane: fit$opt has no record; counted by the unsettled warning),
# and the warnings.
#   Rscript dev/quadseam-acc.R <lib|base> [S]
# log dev/quadseam-log/acc-<tag>.txt
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
S <- if (length(a) > 1) as.integer(a[2]) else 10L
suppressMessages({library(frmtmb); library(lme4)})
cat("lib: frmtmb", format(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
gen <- function(fam, s) {
  set.seed(s)
  if (fam == "bern") {
    d <- data.frame(g = factor(rep(1:100, each = 4)), x = rnorm(400))
    d$y <- rbinom(400, 1, plogis(-0.5 + 0.7 * d$x + rnorm(100, 0, 1.2)[d$g]))
    return(list(d = d, ff = y ~ x + (1 | g), fam = bernoulli(),
                gf = binomial()))
  }
  d <- data.frame(g = factor(rep(1:60, each = 6)), x = rnorm(360))
  eta <- 0.2 + 0.5 * d$x + rnorm(60, 0, 0.3)[d$g]
  if (fam == "bin") {
    d$k <- rbinom(360, 10, plogis(eta)); d$nk <- 10 - d$k
    list(d = d, ff = cbind(k, nk) ~ x + (1 | g), fam = binomial(),
         gf = binomial())
  } else {
    d$y <- rpois(360, exp(eta))
    list(d = d, ff = y ~ x + (1 | g), fam = poisson(), gf = poisson())
  }
}
for (fam in c("bern", "bin", "pois")) {
  dl <- numeric(0)
  for (s in seq_len(S)) {
    g <- gen(fam, s)
    w <- character()
    f <- withCallingHandlers(
      frm(g$ff, family = g$fam, data = g$d, quadrature = TRUE),
      warning = function(x) {
        w <<- c(w, class(x)[1]); invokeRestart("muffleWarning")
      },
      message = function(x) invokeRestart("muffleMessage"))
    l25 <- suppressMessages(suppressWarnings(
      glmer(g$ff, data = g$d, family = g$gf, nAGQ = 25)))
    d1 <- as.numeric(logLik(f)) - as.numeric(logLik(l25))
    dl <- c(dl, d1)
    rs <- f$opt[["quad_settle"]]
    cat(sprintf(paste0("FIT %-4s seed %2d ll %.10f d_glmer25 %10.3e ",
                       "rounds %s grad %s warn %s\n"),
                fam, s, as.numeric(logLik(f)), d1,
                if (is.null(rs)) "-" else rs[["rounds"]],
                if (is.null(rs)) "-" else format(rs[["grad"]], digits = 3),
                if (length(w)) paste(w, collapse = ",") else "-"))
  }
  cat(sprintf("SUM %-4s n %d max|d| %.3e\n", fam, length(dl), max(abs(dl))))
}
cat("END\n")
