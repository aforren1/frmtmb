# Review of lane quadseam: boundary verdicts under one-block quadrature
# on interior and near-zero group sds with few groups, 50 seeds per
# design. Per seed: glmer(nAGQ = 25) sd and isSingular(); the lane
# build's quadrature fit (boundary message, se_lost, NaN SE of the sd,
# sd, logLik, seconds) and its Laplace fit (boundary message, sd).
#   Rscript quadseam-review-bnd.R <base|lane> <design,...> <seeds>
a <- commandArgs(TRUE)
arm <- a[1]
designs <- strsplit(a[2], ",")[[1]]
seeds <- eval(parse(text = a[3]))
base <- "C:/Users/adf44/source/r/rellib-r20"
lane <- "C:/Users/adf44/source/r/wt-quadseam-lib"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (arm == "lane") c(lane, base, ul) else c(base, ul))
suppressMessages({library(frmtmb); library(lme4)})
cat("lib:", arm, dirname(find.package("frmtmb")), "\n")
gen <- function(des, s) {
  set.seed(s)
  mk <- function(G, m) data.frame(g = factor(rep(seq_len(G), each = m)),
                                  x = rnorm(G * m))
  switch(des,
    bin8 = { d <- mk(8, 5); u <- rnorm(8, 0, 0.15); d$N <- 10
      d$k <- rbinom(nrow(d), 10, plogis(0.2 + 0.5 * d$x + u[d$g]))
      d$nk <- d$N - d$k
      list(d = d, ff = cbind(k, nk) ~ x + (1 | g), fam = binomial(),
           gf = binomial()) },
    bern10 = { d <- mk(10, 20); u <- rnorm(10, 0, 0.2)
      d$y <- rbinom(nrow(d), 1, plogis(-0.3 + 0.5 * d$x + u[d$g]))
      list(d = d, ff = y ~ x + (1 | g), fam = bernoulli(), gf = binomial()) },
    pois10 = { d <- mk(10, 10); u <- rnorm(10, 0, 0.1)
      d$y <- rpois(nrow(d), exp(1 + 0.3 * d$x + u[d$g]))
      list(d = d, ff = y ~ x + (1 | g), fam = poisson(), gf = poisson()) },
    bern15s = { d <- mk(15, 10); u <- rnorm(15, 0, 0.3)
      d$y <- rbinom(nrow(d), 1, plogis(0.5 * d$x + u[d$g]))
      list(d = d, ff = y ~ x + (1 | g), fam = bernoulli(), gf = binomial()) },
    stop(des))
}
run <- function(expr) {
  ms <- character(); w <- character()
  t0 <- proc.time()[[3]]
  f <- tryCatch(withCallingHandlers(expr,
    warning = function(x) {
      w <<- c(w, class(x)[1]); invokeRestart("muffleWarning")
    },
    message = function(x) {
      ms <<- c(ms, class(x)[1]); invokeRestart("muffleMessage")
    }), error = function(e) e)
  list(f = f, w = w, ms = ms, secs = proc.time()[[3]] - t0)
}
for (des in designs) for (s in seeds) {
  g <- gen(des, s)
  m <- suppressMessages(suppressWarnings(
    glmer(g$ff, data = g$d, family = g$gf, nAGQ = 25)))
  q <- run(frm(g$ff, family = g$fam, data = g$d, quadrature = TRUE))
  l <- run(frm(g$ff, family = g$fam, data = g$d, quadrature = FALSE))
  info <- function(r) {
    if (inherits(r$f, "error")) return("ERROR")
    sdr <- tryCatch(frmtmb:::sdr_of(r$f), error = function(e) NULL)
    vc <- tryCatch(VarCorr(r$f)[["g"]]$sd, error = function(e) NULL)
    sprintf("sd %.4g ll %.10f bmsg %d lost %s seNaN %s q25 %.3g w %s",
            exp(r$f$estimates$theta[1]), as.numeric(logLik(r$f)),
            sum(r$ms == "frmtmb_boundary_fit"),
            if (is.null(sdr$se_lost)) "-" else
              paste(sdr$se_lost, collapse = ","),
            if (is.null(vc)) "?" else is.nan(vc[1L, "Est.Error"]),
            if (is.null(vc)) NA else vc[1L, "Q2.5"],
            if (length(r$w)) paste(r$w, collapse = ",") else "-")
  }
  cat(sprintf("ROW %s %d glmer sd %.4g ll %.10f sing %s | QUAD %s secs %.2f | LAP %s\n",
              des, s, getME(m, "theta"), as.numeric(logLik(m)),
              isSingular(m), info(q), q$secs, info(l)))
}
cat("END\n")
