# Reviewer, fix round 1: does prior_spike_verdict() exempt only priors
# that bound the dispersion? For every exempted arm the fit must stop at
# a finite dispersion; every refused arm is checked on rellib-r10 for a
# runaway.  Data of dev/fuzzfix-rev-wu.R (seed 20261021, n = 80,
# sum(w) - n = 19.80) plus a small-excess copy (sum(w) - n = 0.5).
#   Rscript dev/fuzzfix-rev-r1-priors.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
run <- function(name, expr) {
  ws <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) e)
  if (inherits(r, "error")) {
    m <- conditionMessage(r)
    tag <- if (grepl("has no maximum|is improper", m)) "REFUSED" else "ERROR"
    cat(sprintf("ARM %-30s %s | %s\n", name, tag, substr(m, 1, 110)))
  } else {
    conv <- r$opt$convergence %||% NA
    cat(sprintf("ARM %-30s ok logLik %.6g | betad %s | conv %s | warnings %d%s\n",
                name, as.numeric(logLik(r)),
                paste(sprintf("%.4g", r$estimates$betad), collapse = " "),
                conv, length(ws),
                if (length(ws)) paste0(" (", substr(ws[1], 1, 60), ")") else ""))
  }
}
`%||%` <- function(a, b) if (is.null(a)) b else a
set.seed(20261021)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n), xs = rnorm(n),
                w = runif(n, 0.5, 2))
eta <- 0.5 + 0.4 * d$x + rnorm(n, 0, 0.5)
d$y <- eta + rnorm(n, 0, exp(0.3 * d$xs))
d$yg <- rgamma(n, shape = 3, rate = 3 / exp(eta))
d$yb <- plogis(eta + rnorm(n, 0, 0.5))
# small excess: sum(w) = n + 0.5
d$w05 <- d$w / sum(d$w) * (n + 0.5)
cat(sprintf("excess big %.4f small %.4f\n", sum(d$w) - n, sum(d$w05) - n))
G <- gaussian()
P <- function(...) set_prior(...)
fo <- bf(y | weights(w) ~ x + (1 | obs))
cat("--- gaussian, sigma without formula (natural scale), excess 19.8\n")
run("halfnormal(0,1) sigma", frm(fo, family = G, data = d,
                                 prior = P("normal(0, 1)", class = "sigma")))
run("exponential(1) sigma", frm(fo, family = G, data = d,
                                prior = P("exponential(1)", class = "sigma")))
run("gamma(2,1) sigma", frm(fo, family = G, data = d,
                            prior = P("gamma(2, 1)", class = "sigma")))
run("gamma(25,1) sigma (order 24)", frm(fo, family = G, data = d,
                             prior = P("gamma(25, 1)", class = "sigma")))
run("lognormal(0,1) sigma", frm(fo, family = G, data = d,
                                prior = P("lognormal(0, 1)", class = "sigma")))
run("inv_gamma(2,1) sigma", frm(fo, family = G, data = d,
                                prior = P("inv_gamma(2, 1)", class = "sigma")))
run("student_t lb=0.05 sigma", frm(fo, family = G, data = d,
                                   prior = P("student_t(3, 0, 2.5)",
                                             class = "sigma", lb = 0.05)))
cat("--- gaussian, small excess 0.5\n")
fo5 <- bf(y | weights(w05) ~ x + (1 | obs))
run("none, excess 0.5", frm(fo5, family = G, data = d))
run("gamma(2,1) sigma, excess 0.5", frm(fo5, family = G, data = d,
                                prior = P("gamma(2, 1)", class = "sigma")))
run("halfnormal sigma, excess 0.5", frm(fo5, family = G, data = d,
                                prior = P("normal(0, 1)", class = "sigma")))
cat("--- gaussian, sigma ~ 1 + xs (log-scale intercept)\n")
fos <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + xs)
run("normal(0,1) Intercept_sigma", frm(fos, family = G, data = d,
    prior = P("normal(0, 1)", class = "Intercept", dpar = "sigma")))
run("student_t Intercept_sigma", frm(fos, family = G, data = d,
    prior = P("student_t(3, 0, 2.5)", class = "Intercept", dpar = "sigma")))
run("normal b_sigma only", frm(fos, family = G, data = d,
    prior = P("normal(0, 1)", class = "b", dpar = "sigma")))
cat("--- Gamma shape\n")
fg <- bf(yg | weights(w) ~ x + (1 | obs))
GA <- Gamma("log")
run("gamma(0.01,0.01) shape", frm(fg, family = GA, data = d,
                                  prior = P("gamma(0.01, 0.01)", class = "shape")))
run("student_t(3,0,2.5) shape", frm(fg, family = GA, data = d,
                                    prior = P("student_t(3, 0, 2.5)", class = "shape")))
run("exponential(1) shape", frm(fg, family = GA, data = d,
                                prior = P("exponential(1)", class = "shape")))
run("normal(0,10) shape", frm(fg, family = GA, data = d,
                              prior = P("normal(0, 10)", class = "shape")))
cat("--- Beta phi\n")
fb <- bf(yb | weights(w) ~ x + (1 | obs))
run("none phi", frm(fb, family = Beta(), data = d))
run("gamma(0.01,0.01) phi", frm(fb, family = Beta(), data = d,
                                prior = P("gamma(0.01, 0.01)", class = "phi")))
run("student_t(3,0,2.5) phi", frm(fb, family = Beta(), data = d,
                                  prior = P("student_t(3, 0, 2.5)", class = "phi")))
cat("DONE\n")
