source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
set.seed(1)
n <- 200
x <- rnorm(n) * 1e-4 + 5e-4
dg <- data.frame(x = x, y = 1 + 2000 * (x - 5e-4) + rnorm(n))
for (a in list(NULL, FALSE, TRUE)) {
  f <- suppressWarnings(frm(bf(y ~ x), family = gaussian(), data = dg,
    prior = set_prior("", class = "Intercept", ub = 0.5),
    control = frmtmb_control(autoscale = a)))
  cat("autoscale", format(a), "\n"); print(f$opt$par, digits = 12)
  print(f$obj$gr(f$opt$par)); print(f$par_units); cat(f$opt$message, "\n")
}
