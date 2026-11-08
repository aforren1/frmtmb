# Review of lane nlstart: the fit's state after the sign probe.
# Usage: REV_ARM=after|<mutant lib> Rscript dev/review-nlstart-state.R
# (1) Fits the probe cannot improve (no winner): every field against the
#     probe-off fit, identical() or the largest difference.
# (2) A winner: is every field the winner's?
source("dev/review-nlstart-common.R")
same <- function(a, b) {
  if (identical(a, b)) return("identical")
  if (is.numeric(a) && is.numeric(b) && length(a) == length(b)) {
    return(sprintf("DIFF max %.3g", max(abs(a - b))))
  }
  "DIFF"
}
fields <- function(f) {
  env <- f$obj$env
  lf <- env$lfixed()
  sdr <- frmtmb:::sdr_of(f)
  list(opt_par = f$opt$par, opt_obj = f$opt$objective,
       conv = f$opt$convergence, msg = f$opt$message,
       last_par = env$last.par, last_par_best = env$last.par.best,
       value_best = env$value.best,
       sdr_par = sdr$par.fixed, sdr_cov = sdr$cov.fixed,
       est = unlist(f$estimates), start = f$start, units = f$par_units,
       vcov = suppressWarnings(vcov(f)),
       fitted = suppressWarnings(fitted(f)),
       ranef = tryCatch(unlist(ranef(f)), error = function(e) NULL))
}
compare <- function(label, fm, d, fam = gaussian()) {
  on <- tt(frm(fm, family = fam, data = d))
  off <- withr::with_options(list(frmtmb.nl_sign_probe = FALSE),
                             tt(frm(fm, family = fam, data = d)))
  cat("==", label, "on:", llv(on), "off:", llv(off), "multistart:",
      any(attr(on, "mc") == "frmtmb_nl_multistart"), "\n")
  a <- fields(on); b <- fields(off)
  for (nm in names(a)) cat(sprintf("  %-14s %s\n", nm, same(a[[nm]], b[[nm]])))
  # what a reader that evaluates the objective at the best point sees
  pa <- on$obj$env$last.par.best
  cat("  on: last.par.best[lfixed] == opt$par:",
      identical(unname(pa[on$obj$env$lfixed()]), unname(on$opt$par)), "\n")
  invisible(list(on = on, off = off))
}
set.seed(7302)
n <- 300
d <- data.frame(x = round(rnorm(n), 3))
d$y <- 2 * exp(0.5 * d$x) + rnorm(n, 0, 0.3)
compare("exp, no RE (test's data, seed 7302)",
        bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), d)
set.seed(7311)
d <- data.frame(x = round(rnorm(n), 3), g = factor(sample(1:12, n, TRUE)))
d$y <- (2 + rnorm(12, 0, 0.5)[d$g]) * exp(0.5 * d$x) + rnorm(n, 0, 0.3)
compare("exp, b1 ~ 1 + (1 | g), seed 7311",
        bf(y ~ b1 * exp(b2 * x), b1 ~ 1 + (1 | g), b2 ~ 1, nl = TRUE), d)
set.seed(7312)
d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
d$y <- 2 * exp(0.5 * d$x) + rnorm(n, 0, exp(-1 + 0.4 * d$z))
compare("exp, sigma ~ z, seed 7312",
        bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, sigma ~ z, nl = TRUE), d)

cat("\n== winner: sat body, data of test-nl-start.R seed 101\n")
set.seed(101)
d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
d$y <- 1.5 * exp(0.6 * d$x) + 0.3 * d$z + rnorm(n, 0, 0.5)
fm <- bf(y ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE)
w <- tt(frm(fm, data = d))
cat("logLik", llv(w), "conv", w$opt$convergence, "messages:",
    paste(attr(w, "mc"), collapse = ","), "warnings:", length(attr(w, "w")),
    "\n")
for (m in attr(w, "m")) cat("  MSG:", m, "\n")
env <- w$obj$env
cat("last.par.best[lfixed] == opt$par:",
    identical(unname(env$last.par.best[env$lfixed()]), unname(w$opt$par)), "\n")
cat("-opt$objective vs logLik:", -w$opt$objective - llv(w), "\n")
sdr <- frmtmb:::sdr_of(w)
cat("sdr par.fixed vs opt$par:", same(unname(sdr$par.fixed), unname(w$opt$par)), "\n")
cat("fit$start is NULL:", is.null(w$start), "\n")
if (!is.null(w$start)) {
  cat("fit$start beta:", format(w$start$beta, digits = 8), "\n")
}
# a fit started where the probe's winner was refitted from
s <- tt(frm(fm, data = d, start = w$start))
cat("frm(start = fit$start):", llv(s), "same opt$par:",
    same(unname(s$opt$par), unname(w$opt$par)), "\n")
# frm_allfit(start_from_mle = FALSE) reads fit$start
af <- tryCatch(frm_allfit(w, start_from_mle = FALSE),
               error = function(e) conditionMessage(e))
if (is.character(af)) cat("allfit error:", af, "\n") else {
  print(af)
}
# and on the probe-off fit, from the same data (default start NULL)
woff <- withr::with_options(list(frmtmb.nl_sign_probe = FALSE),
                            tt(frm(fm, data = d)))
cat("probe-off fit$start is NULL:", is.null(woff$start), "\n")
print(tryCatch(frm_allfit(woff, start_from_mle = FALSE),
               error = function(e) conditionMessage(e)))
up <- tt(update(w))
cat("update(fit):", llv(up), "multistart:",
    any(attr(up, "mc") == "frmtmb_nl_multistart"), "\n")
