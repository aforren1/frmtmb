# Reviewer of lane fixes, claim 2: what stats::D() returns for the
# calls a nonlinear body uses, and where "identical derivatives" does
# not mean "a function of the sum".
tryD <- function(e, v) tryCatch(deparse1(stats::D(e, v)),
                                error = function(err) paste("ERROR:", conditionMessage(err)))
for (e in list(quote(pnorm(pt, m2, s)), quote(pnorm(pt, m1, s) + pnorm(pt, m2, s)),
               quote(dnorm(x, a, s) + dnorm(x, b, s)), quote(log(a, b)),
               quote(psigamma(a, b)), quote(besselJ(a, b)),
               quote(a + b), quote(b + a), quote((a + b) * x), quote(a + 2 * b),
               quote(a - b), quote(exp(a + b)), quote(inv_logit(a + b)),
               quote(ifelse(x > 0, a, b)), quote(a * 0 + b * 0 + c),
               quote(pmax(a, b)), quote(round(a) + round(b)), quote(a %% 1 + b),
               quote(sign(a) + sign(b)), quote(abs(a) + abs(b)))) {
  cat(sprintf("%-42s D_a: %-30s D_b: %s   D_m1: %s  D_m2: %s\n", deparse1(e),
              tryD(e, "a"), tryD(e, "b"), tryD(e, "m1"), tryD(e, "m2")))
}
