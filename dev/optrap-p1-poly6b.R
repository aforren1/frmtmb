# Lane optrap, punch round 1 (m10): why the degree-6 raw polynomial at
# scale 1 stops 3.33 short at code 0 with a unit of 0.25. Prints the
# stop, then restarts nlminb from it with the same scale, with scale 1,
# and evaluates the gradient there in natural units.
#   Rscript dev/optrap-p1-poly6b.R <lib> [round0]
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
if (identical(a[2], "round0")) source("dev/optrap-round0-unit.R")
deg <- 6
set.seed(deg)
d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
fo <- as.formula(paste("y ~", paste0("x", 1:deg, collapse = " + ")))
l <- lm(fo, data = d)
f <- suppressWarnings(frm(fo, data = d))
o <- f$opt
cat("units", deparse(f$par_units), "\n")
cat("stop: code", o$convergence, o$message, "iter", o$iterations,
    "evals", deparse(o$evaluations), "| logLik - lm",
    -o$objective - as.numeric(logLik(l)), "\n")
u <- f$par_units %||% rep(1, length(o$par))
g <- f$obj$gr(o$par)
cat("gradient x unit, max", max(abs(g * u)), " raw max", max(abs(g)), "\n")
for (sc in list(same = 1 / u, one = rep(1, length(u)))) {
  r <- nlminb(o$par, f$obj$fn, f$obj$gr, scale = sc,
              control = list(eval.max = 5000, iter.max = 5000))
  cat("restart: code", r$convergence, r$message, "iter", r$iterations,
      "| logLik - lm", -r$objective - as.numeric(logLik(l)), "\n")
}
cat("autoscale on this fit:", deparse(f$frame$autoscale %||% "none"), "\n")
