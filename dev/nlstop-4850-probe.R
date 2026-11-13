# Lane nlstop, item 4850: what the resolved prior entries hold for y ~ x
# with normal(0, 1000) on b at y about 1e6, and where the optimizer goes
# from two restarts with b at its prior location: the rest at the tape's
# start, and the rest at the fit.
source("dev/nlstop-common.R")
ns <- asNamespace("frmtmb")
set.seed(6201)
n <- 200
d <- data.frame(x = rnorm(n))
d$y <- 0.5 * d$x + rnorm(n)
for (s in c(1e-2, 1e-4, 1e-8, 1e-20)) for (psd in 1:3) invisible(NULL)
s <- 1e6
d$yl <- s * (1 + 0.05 * d$x) + rnorm(n, 0, s * 1e-3)
pr_l <- c(set_prior("normal(0, 1000)", class = "b"),
          set_prior(sprintf("normal(%g, %g)", s, s), class = "Intercept"))
f <- frm(yl ~ x, data = d, prior = pr_l)
ri <- ns$resolve_prior_input(list(frame = f$frame, spec = f$spec), pr_l)
str(ri$entries, max.level = 3)
obj <- f$obj
cat("tape start par:", format(obj$par, digits = 8), "\n")
cat("fit par:", format(f$opt$par, digits = 8), " obj", f$opt$objective, "\n")
ctl <- f$control
bounds <- ns$fit_outer_box(f)
for (lab in c("start", "fit")) {
  p <- if (lab == "start") obj$par else f$opt$par
  p[2] <- 0
  o <- ns$optimize_obj(obj, ctl, bounds, f$par_units, start_par = p)
  cat(lab, ": objective", format(o$objective, digits = 10), "par",
      format(o$par, digits = 8), "code", o$convergence, "\n")
}
cat("par_units:", format(f$par_units, digits = 4), "\n")
