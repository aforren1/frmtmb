# Review of lane quadseam2, attack 2: which one-block quadrature fits
# does quad_sd_up_escape() touch, and what does it cost? The design of
# dev/quadseam2-break.R (15 groups of 4, true group sd 0 and 0.05, six
# families), where most fits end with the sd near zero and the escape's
# probe runs. Per fit: logLik in hex, the sd, sd_escape, fit and vcov
# seconds, warnings and SE verdict.
#   Rscript quadseam2-review-touch.R <base|lib>
a <- commandArgs(TRUE)
if (a[1] == "lib") a[1] <- "C:/Users/adf44/source/r/wt-quadseam2-lib"
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
`%||%` <- function(a, b) if (is.null(a)) b else a
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
for (fn in names(fams)) for (sdg in c(0, 0.05)) for (s in 1:10) {
  d <- gen(fn, s, sdg)
  w <- character()
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(withCallingHandlers(
    frm(y ~ x + (1 | g), family = fams[[fn]], data = d, quadrature = TRUE),
    warning = function(x) {
      w <<- c(w, class(x)[1]); invokeRestart("muffleWarning")
    }, message = function(x) {
      w <<- c(w, paste0("m:", class(x)[1])); invokeRestart("muffleMessage")
    }), error = function(e) e)
  t1 <- proc.time()[["elapsed"]]
  if (inherits(f, "error")) {
    cat(sprintf("ROW\t%s\t%.2f\t%d\tERROR\n", fn, sdg, s))
    next
  }
  se <- tryCatch(withCallingHandlers(sqrt(diag(vcov(f))),
    warning = function(x) {
      w <<- c(w, paste0("v:", class(x)[1])); invokeRestart("muffleWarning")
    }, message = function(x) {
      w <<- c(w, paste0("vm:", class(x)[1])); invokeRestart("muffleMessage")
    }), error = function(e) NA)
  t2 <- proc.time()[["elapsed"]]
  cat(sprintf("ROW\t%s\t%.2f\t%d\t%s\t%.12f\t%.3g\t%s\t%.2f\t%.2f\t%s\t%s\n",
              fn, sdg, s, sprintf("%a", as.numeric(logLik(f))),
              as.numeric(logLik(f)), exp(f$estimates$theta[1]),
              format(f$opt$sd_escape[["gain"]] %||% NA, digits = 3),
              t1 - t0, t2 - t1, paste(sort(unique(w)), collapse = ","),
              paste(names(se)[!is.finite(se)], collapse = ",")))
}
cat("END\n")
