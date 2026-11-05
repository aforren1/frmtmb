# exp(0.8 x), seeds 5701, 5702: does restarting from a false-convergence stop climb
# the ridge to the closed form's maximum? nlminb from the stop, ten times, against a
# Newton step with the closed form's own Hessian for comparison of the path.
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
for (seed in c(5701, 5702)) {
  set.seed(seed)
  n <- 100
  x <- sort(runif(n, 0, 10))
  e0 <- rnorm(n)
  y <- exp(0.8 * x) + 0.3 * e0
  d <- data.frame(x = x, y = y)
  fit <- suppressWarnings(frm(bf(y ~ gp(x)), data = d, family = gaussian()))
  obj <- fit$obj
  p <- fit$opt$par
  cat("seed", seed, "stop", -fit$opt$objective, fit$opt$message, "\n")
  un <- fit$par_units %||% rep(1, length(p))
  for (k in 1:6) {
    o <- nlminb(p, obj$fn, obj$gr, control = list(iter.max = 300, eval.max = 400))
    cat(sprintf(" restart %d: %s ll %.5f th %s\n", k, o$message, -o$objective,
                paste(sprintf("%.5f", o$par), collapse = " ")))
    if (o$objective >= min(obj$env$value.best, Inf) + 1e-6 && k > 1) NULL
    p <- o$par
  }
  # from theta on a grid of length scales, the same sd and sigma
  for (rho in c(3.5, 3.8, 4.0, 4.1, 4.14, 4.2, 4.4, 4.8)) {
    q <- fit$opt$par
    q[which(names(q) == "theta")[2]] <- log(rho)
    o <- tryCatch(nlminb(q, obj$fn, obj$gr, control = list(iter.max = 300, eval.max = 400)), error = function(e) NULL)
    if (!is.null(o)) cat(sprintf(" start rho %.2f: %s ll %.5f rho %.4f sd %.0f\n", rho, o$message, -o$objective, exp(o$par[names(o$par) == "theta"][2]), exp(o$par[names(o$par) == "theta"][1])))
  }
}
