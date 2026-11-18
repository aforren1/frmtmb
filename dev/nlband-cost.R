# Lane nlband: what the hole start costs a fit. One process; the arms are
# "on" (the lane's escape), "off" (nl_pole_holes() returning NULL, which
# is 0.83.0's escape) and "off2" (the same as off: the control, which
# must read 1 against off). Each block fits its rows twice; the arms
# rotate within a round and the minimum of 3 rounds is kept. Then the
# load-independent count: the objective and gradient evaluations of
# every optimize_obj() run of each fit, on and off.
# NLBAND_ARM=after Rscript dev/nlband-cost.R
source("dev/nlband-common.R")
holes_fn <- ns$nl_pole_holes
set_holes <- function(on) {
  f <- if (on) holes_fn else function(w) NULL
  unlockBinding("nl_pole_holes", ns)
  assign("nl_pole_holes", f, envir = ns)
  lockBinding("nl_pole_holes", ns)
}
q <- function(expr) suppressWarnings(suppressMessages(expr))
worse <- function(s, n = 400) {
  d <- band_data(s, n = n)
  d$y <- 1 / (d$x + 0.5) + rnorm(nrow(d), 0, 0.2)
  d
}
nohole <- function(s, n = 400) {
  set.seed(s)
  d <- data.frame(x = runif(n, 0, 4))
  d$y <- 1 / (d$x + 0.5) + rnorm(n, 0, 0.2)
  d
}
jitter <- function(s) {
  set.seed(s)
  x <- rep(1:60, each = 3) + runif(180, -0.01, 0.01)
  data.frame(x = x, y = 5 / (x + 0.5) + rnorm(180, 0, 0.2))
}
sets <- list(
  "hole wins (band data, n 400)" = lapply(1:3, band_data),
  "hole loses (truth b 0.5, n 400)" = lapply(1:3, worse),
  "no hole (n 400)" = lapply(1:3, nohole),
  "hole loses (truth b 0.5, n 20000)" = lapply(1:2, worse, n = 20000),
  "two holes lose (jittered replicates, review)" = lapply(1:3, jitter))
arms <- c("on", "off", "off2")
for (nm in names(sets)) {
  rows <- sets[[nm]]
  tmin <- setNames(rep(Inf, 3), arms)
  for (round in 1:3) {
    for (a in arms[c(round:3, seq_len(round - 1))]) {
      set_holes(a == "on")
      t0 <- proc.time()[["elapsed"]]
      for (rep in 1:2) for (d in rows) q(frm(band_form, data = d))
      tmin[a] <- min(tmin[a], proc.time()[["elapsed"]] - t0)
    }
  }
  set_holes(TRUE)
  cat(sprintf("COST %s: min s on %.2f off %.2f off2 %.2f | on/off %.3f control %.3f\n",
              nm, tmin["on"], tmin["off"], tmin["off2"],
              tmin["on"] / tmin["off"], tmin["off2"] / tmin["off"]))
}
n_ev <- 0
suppressMessages(trace("optimize_obj",
                       exit = quote(n_ev <<- n_ev + returnValue()$evals),
                       print = FALSE, where = ns))
for (nm in names(sets)) for (i in seq_along(sets[[nm]])) {
  ev <- c(on = NA, off = NA)
  for (a in names(ev)) {
    set_holes(a == "on")
    n_ev <- 0
    q(frm(band_form, data = sets[[nm]][[i]]))
    ev[a] <- n_ev
  }
  set_holes(TRUE)
  cat(sprintf("COUNT %s row %d: evaluations on %d off %d ratio %.3f\n",
              nm, i, as.integer(ev[["on"]]), as.integer(ev[["off"]]),
              ev[["on"]] / ev[["off"]]))
}
cat("DONE\n")
