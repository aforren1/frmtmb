# Reviewer of lane ordmix, re-check: does the new non-finite-gradient
# warning of check_convergence() fire on a CORRECT optimum where a
# parameter sits at the edge of its domain? Each case is a fit whose
# maximum is at or toward a boundary: a box bound that holds, a
# variance component collapsing to 0, a zero-inflation with no zeros,
# overdispersion with none in the data, student-t on normal data, a
# kinked (asymmetric Laplace) likelihood, a held dpar, profile = TRUE,
# REML, a hurdle with no zeros. Seed 20261006. Every warning recorded.
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261006)
n <- 300
d <- data.frame(x = rnorm(n), g = factor(sample(1:15, n, TRUE)))
d$y <- 1 + 0.5 * d$x + rnorm(n)
d$cnt <- rpois(n, exp(0.3 + 0.2 * d$x))
d$cnt1 <- d$cnt + 1L
d$ybin <- rbinom(n, 1, plogis(0.3 * d$x))
d$o <- 1L + (0.8 * d$x + rlogis(n) > c(-1)) + (0.8 * d$x + rlogis(n) > 1)
run <- function(lab, expr) {
  w <- character(0)
  f <- tryCatch(withCallingHandlers(expr, warning = function(e) {
    w <<- c(w, conditionMessage(e))
    invokeRestart("muffleWarning")
  }), error = function(e) conditionMessage(e))
  if (is.character(f)) {
    cat(sprintf("%-28s ERROR %s\n", lab, substr(f, 1, 120)))
    return(invisible())
  }
  g <- tryCatch(f$obj$gr(f$opt$par), error = function(e) NA)
  cat(sprintf("%-28s gradfinite=%s maxgrad=%.3g nonfinite_warning=%s other_warnings=%d\n",
              lab, all(is.finite(g)), max(abs(g)),
              any(grepl("gradient at the reported optimum is not finite", w,
                        fixed = TRUE)),
              sum(!grepl("gradient at the reported optimum is not finite", w,
                         fixed = TRUE))))
}
run("box bound holds (lb = 1)",
    frm(bf(y ~ x), data = d,
        prior = set_prior("", class = "b", coef = "x", lb = 1)))
run("RE variance -> 0",
    frm(bf(y ~ x + (1 | g)), data = d))
run("zi poisson, no zeros",
    frm(bf(cnt1 ~ x), family = zero_inflated_poisson(), data = d))
run("negbinomial, poisson data",
    frm(bf(cnt ~ x), family = negbinomial(), data = d))
run("student, normal data",
    frm(bf(y ~ x), family = student(), data = d))
run("asym_laplace (kink)",
    frm(bf(y ~ x), family = asym_laplace(), data = d))
run("held sigma",
    frm(bf(y ~ x, sigma = 1), data = d))
run("profile = TRUE",
    frm(bf(y ~ x + (1 | g)), data = d,
        control = frmtmb_control(profile = TRUE)))
run("REML",
    frm(bf(y ~ x + (1 | g)), data = d, REML = TRUE))
run("hurdle_poisson, no zeros",
    frm(bf(cnt1 ~ x), family = hurdle_poisson(), data = d))
run("bernoulli, separated x",
    frm(bf(sep ~ x), family = bernoulli(),
        data = transform(d, sep = as.integer(x > 0))))
run("cumulative thres(4) above",
    frm(bf(o | thres(3) ~ x), family = cumulative(), data = d))
run("cumulative, an empty middle",
    frm(bf(oe ~ x), family = cumulative(),
        data = transform(d, oe = ifelse(o == 2, 1L, o))))
