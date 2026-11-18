# Lane nlband: why does the escape not reach its loop on the limit
# battery's mm rows? Rscript dev/nlband-dbg1.R
source("dev/nlband-common.R")
set.seed(1)
n <- 50
d <- data.frame(x = runif(n, 0, 10))
d$y <- 10 * d$x / (1 + d$x) + rnorm(n, 0, 0.5)
f <- frm(bf(y ~ Vm * x / (K + x), Vm + K ~ 1, nl = TRUE), data = d)
ctl <- f$control
cat("nl_starts", format(ctl$nl_starts), "quad", format(f$quadrature),
    "imp", is.null(f$importance), "whiten", is.null(f$beta_whiten),
    "profile", format(ctl$profile), "optimizer", format(ctl$optimizer),
    "stops", ns$optctrl_stops_earlier(ctl$optCtrl), "\n")
lp <- Filter(function(q) !is.null(q$nl_body), f$frame$linpreds)[[1L]]
wl <- ns$nl_pole_walls(f$frame, lp, f$estimates, "beta", 2L,
                       ns$nl_denoms(lp$nl_body))
str(wl)
print(ns$nl_pole_escape(function(...) stop("no refit"), f$frame, f))
