# Item 5: NEWS says the L1 constant on sds (-0.667, -0.848, -1.97) is because
# "brms normalizes the half density". Measure it: with and without the prior,
# per form. The half density would add log 2 per standard deviation, 0.693.
# Seed 10400.
source("dev/review-refuse-lib.R")
set.seed(10400)
n <- 200
d <- data.frame(x = runif(n), z = runif(n), f = factor(sample(c("a", "b", "c"), n, TRUE)))
d$y <- sin(6 * d$x) + d$z^2 + rnorm(n, 0, 0.4)
p1 <- "set_prior('student_t(3, 0, 1)', class = 'sds')"
one <- function(lab, form, nsd) {
  a <- tr(meas(list(bf = form, data = d)))
  b <- tr(meas(list(bf = form, data = d, prior = p1)))
  ca <- if (is_err(a)) NA else (a$l1_const %||% NA)
  cb <- if (is_err(b)) NA else (b$l1_const %||% NA)
  pr("%-22s flat: const %s pts[%s] err[%s]; prior: const %s pts[%s]; prior - flat %.4f, %d sd x log 2 = %.4f",
     lab, format(ca, digits = 5),
     if (is_err(a)) msg(a, 60) else paste(format(unlist(a$l1_extra_const), digits = 3), collapse = ","),
     if (is_err(a)) "" else substr(a$l1_error %||% "", 1, 40), format(cb, digits = 5),
     if (is_err(b)) msg(b, 60) else paste(format(unlist(b$l1_extra_const), digits = 3), collapse = ","),
     cb - ca, nsd, nsd * log(2))
}
one("s(x)", "bf(y ~ s(x))", 1L)
one("t2(x, z)", "bf(y ~ t2(x, z))", 3L)
one("s(x, by = f)", "bf(y ~ f + s(x, by = f))", 3L)
one("s(x) + s(z, k = 5)", "bf(y ~ s(x) + s(z, k = 5))", 2L)
