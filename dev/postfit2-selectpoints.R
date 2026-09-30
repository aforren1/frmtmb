# Why brms and frmtmb keep different points under select_points > 0:
# which variables brms's make_point_frame() measures the distance on.
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-postfit2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(brms)
  library(frmtmb)
})
set.seed(1)
n <- 150
d <- data.frame(x = runif(n), z = runif(n),
                f = factor(sample(c("a", "b"), n, TRUE)))
d$y <- sin(2 * pi * d$x) + d$z * (d$f == "b") + rnorm(n, 0, 0.3)
b3 <- suppressMessages(suppressWarnings(
  brm(y ~ f + x + z, data = d, algorithm = "fixed_param", chains = 1,
      iter = 1, warmup = 0, refresh = 0, seed = 1, silent = 2,
      init = list(list(b = array(c(0.1, 0.2, 0.3), 3), Intercept = 0.2,
                       sigma = 1)))))
cond <- brms:::prepare_conditions(b3, effects = list("f"))
print(cond)
print(names(model.frame(b3)))
ce <- conditional_effects(b3, "f", select_points = 0.1)
pts <- attr(ce[[1]], "points")
print(pts)
unit <- function(v, at) abs((v - min(v)) / diff(range(v)) -
                              (at - min(v)) / diff(range(v)))
k_xz <- unit(d$x, mean(d$x)) <= 0.1 & unit(d$z, mean(d$z)) <= 0.1
cat("rows within 0.1 on x and z:", sum(k_xz), "\n")
for (yv in unique(cond$y)) {
  k_y <- k_xz & unit(d$y, yv) <= 0.1
  cat("... and on y at", yv, ":", sum(k_y), "\n")
}
fce <- frmtmb::conditional_effects(frm(bf(y ~ f + x + z), data = d), "f",
                                   select_points = 0.1)
print(attr(fce[[1]], "points"))
print(fce[[1]][1, ])
