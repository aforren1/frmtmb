source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-helpers.R")
# B1: make_start()'s announced warning at a bounded link. The reported
# fit reuses the pre-fit's template and never calls make_start(), so an
# engaged default cannot announce it. Seed 507.
set.seed(507)
n <- 200
d <- data.frame(x = rnorm(n) * 1e-4, y0 = 0L, xs = rnorm(n))
run <- function(lab, ...) {
  w <- character(0)
  f <- withCallingHandlers(tryCatch(frm(...), error = function(e) e),
         warning = function(x) { w <<- c(w, conditionMessage(x))
                                 invokeRestart("muffleWarning") })
  cat(sprintf("%-28s %s | %d warning(s)\n", lab,
              if (inherits(f, "error")) conditionMessage(f) else
                sprintf("logLik %.4g code %d", as.numeric(logLik(f)),
                        f$opt$convergence), length(w)))
  for (x in w) cat("     |", substr(x, 1, 140), "\n")
}
run("poisson y=0, x 1e-4, default", y0 ~ x, family = poisson(), data = d)
run("poisson y=0, x 1e-4, FALSE", y0 ~ x, family = poisson(), data = d,
    control = frmtmb_control(autoscale = FALSE))
run("poisson y=0, x sd 1, default", y0 ~ xs, family = poisson(), data = d)
