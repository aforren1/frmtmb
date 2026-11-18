# Review of lane nbrobust, check 4: the negative-binomial CDF that cens()
# and trunc() read (nbinom_cdf_log(), not changed by the lane) at a large
# shape. (a) Its value and derivative in the log shape on the tape against
# a 300-bit Rmpfr reference (sum of the textbook masses; derivative by
# central difference of the mpfr value). (b) Censored and truncated fits
# of item 4903's data at the poisson limit.
#   Rscript dev/reviews/nbrobust-review-cens.R <base|lib> <a|b>
args <- commandArgs(TRUE)
arm <- args[1]
part <- args[2]
base <- c("C:/Users/adf44/source/r/rellib-r22",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lib") c("C:/Users/adf44/source/r/wt-nbrobust-lib", base)
          else base)
suppressPackageStartupMessages({
  library(frmtmb)
  library(Rmpfr)
})
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
cat("La:", La_library(), " CORETYPE:", Sys.getenv("OPENBLAS_CORETYPE"), "\n")
eps <- .Machine$double.eps
P <- 300
m <- function(x) mpfr(x, P)
ref_cdf <- function(q, a, b) {
  mu <- exp(a)
  n <- exp(b)
  k <- 0:q
  lm <- lgamma(k + n) - lgamma(n) - lgamma(m(k + 1)) +
    n * log(n / (n + mu)) + k * log(mu / (n + mu))
  F <- sum(exp(lm))
  c(lF = log(F), lS = log(1 - F))
}
if (part == "a") {
  ncl <- frmtmb:::nbinom_cdf_log
  cat("q lmu lsh | lF err, dlF/dlsh tape, ref, err | lS err, dlS/dlsh tape,",
      "ref, err\n")
  worst <- c(F = 0, S = 0)
  for (q in c(0, 1, 3, 8, 20)) for (lmu in c(0, 1.5)) {
    for (lsh in c(2, 10, 15, 18, 20, 24, 28, 30, 35, 40)) {
      tF <- RTMB::MakeTape(function(p) ncl(q, lmu, p)$lF, lsh)
      tS <- RTMB::MakeTape(function(p) ncl(q, lmu, p)$lS, lsh)
      a <- m(lmu)
      b <- m(lsh)
      h <- m(2)^-90
      r0 <- ref_cdf(q, a, b)
      rd <- (ref_cdf(q, a, b + h) - ref_cdf(q, a, b - h)) / (2 * h)
      gF <- as.numeric(tF$jacobian(lsh))
      gS <- as.numeric(tS$jacobian(lsh))
      eF <- abs(gF - as.numeric(rd[1]))
      eS <- abs(gS - as.numeric(rd[2]))
      worst <- pmax(worst, c(eF, eS))
      cat(sprintf(paste("%2d %3g %3g | %9.2g %11.4g %11.4g %9.2g |",
                        "%9.2g %11.4g %11.4g %9.2g\n"),
                  q, lmu, lsh, abs(tF(lsh) - as.numeric(r0[1])), gF,
                  as.numeric(rd[1]), eF, abs(tS(lsh) - as.numeric(r0[2])),
                  gS, as.numeric(rd[2]), eS))
    }
  }
  cat("worst derivative error: lF", worst[1], " lS", worst[2], "\n")
} else {
  catch <- function(expr) {
    w <- character(0)
    mc <- character(0)
    v <- withCallingHandlers(expr, warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }, message = function(x) {
      mc <<- c(mc, class(x)[1])
      invokeRestart("muffleMessage")
    })
    list(v = v, w = w, mc = mc)
  }
  set.seed(2)
  d0 <- data.frame(x = stats::rnorm(200))
  d0$y <- stats::rpois(200, exp(1 + 0.4 * d0$x))
  d0$cc <- ifelse(d0$y >= 5, "right", "none")
  d0$yc <- pmin(d0$y, 5)
  d0$yt <- d0$y
  for (k in c(0, 4, 8, 12, 16)) {
    d <- d0
    d$x <- d$x * (1 + 2e-15 * k)
    for (form in c("cens", "trunc", "cens30")) {
      r <- catch(switch(form,
        cens = frm(yc | cens(cc) ~ x, family = negbinomial(), data = d),
        cens30 = frm(yc | cens(cc) ~ x, family = negbinomial(), data = d,
                     start = list(betad = 30)),
        trunc = frm(yt | trunc(ub = 12) ~ x, family = negbinomial(),
                    data = d)))
      f <- r$v
      g <- tryCatch(f$obj$gr(f$opt$par), error = function(e) NA)
      cat(sprintf(paste("k %2d %-6s code %d lsh %8.4f logLik %.12f",
                        "gr(shape) %10.3g maxgr %9.3g\n"),
                  k, form, f$opt$convergence, f$opt$par[["betad"]],
                  as.numeric(logLik(f)),
                  g[names(f$opt$par) == "betad"][1], max(abs(g))))
      if (length(r$mc)) cat("   messages:", paste(unique(r$mc), collapse = ", "), "\n")
      for (w in unique(r$w)) cat("   WARNING:", substr(gsub("\n", " ", w), 1, 150), "\n")
    }
  }
}
