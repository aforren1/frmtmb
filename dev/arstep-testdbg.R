# Lane arstep: where the test's formula reference and the tape part, on the
# cases of test-gp-ar-start.R's first test. Args: arm
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/arstep-lib.R")
rev_libs(args[1])
ex <- parse("C:/Users/adf44/source/r/frmtmb-wt-arstep/tests/testthat/test-gp-ar-start.R")
# the helpers only: every top-level assignment
for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e, globalenv())
cases <- list(
  p1 = list(d = ga_step(40L, 2L), f = y ~ gp(x) + ar(time = t), g = FALSE),
  p2 = list(d = ga_step(40L, 1L, p = 2L), f = y ~ gp(x) + ar(time = t, p = 2), g = FALSE),
  gr = list(d = ga_step(48L, 1L, groups = 3L), f = y ~ gp(x) + ar(time = t, gr = g), g = TRUE))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fit <- frm(cs$f, data = cs$d, family = gaussian())
  p0 <- fit$opt$par
  for (k in 0:2) {
    p <- p0 + k * 0.2 * cos(seq_along(p0) + k)
    v <- fit$obj$fn(p)
    nat <- ga_natural(fit, p)
    ref <- ga_ref_nll(cs$d$y, cs$d$x, if (cs$g) cs$d$g else rep(1L, nrow(cs$d)),
                      nat$beta, nat$sigma, nat$sdgp, nat$lscale, nat$ar)
    # the same with a nugget of zero and of 1e-12 * sdgp^2 (absolute)
    cat(sprintf("TD %s k %d tape %.10f ref %.10f diff %.3e | sdgp %.4g lscale %.4g sigma %.4g ar %s\n", nm, k, v, ref, v - ref,
                nat$sdgp, nat$lscale, nat$sigma, paste(signif(nat$ar, 4), collapse = ",")))
  }
}
cat("DONE\n")
