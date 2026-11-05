# The review's M6 row: gp(x, k = 30) on sin(8 x) at noise 1e-3 reads sigma 720
# times high with conv 0. A basis of 30 can carry 7.6 cycles, so is this the
# approximation or the start? The same fit from the theta of k = 60, and the
# frequency the basis reaches. Seed 8102's construction.
#   Rscript dev/gpgap-kstart.R
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(8102)
for (w in c(2, 5, 8, 12)) {
  n2 <- 200
  x <- round(runif(n2, 0, 6), 3)
  e <- rnorm(n2)
  for (s in c(1e-2, 1e-3)) {
    d <- data.frame(x = x, y = sin(w * x) + s * e)
    if (w != 8 || s != 1e-3) next
    fits <- list()
    for (k in c(15, 30, 60, 100)) {
      f <- suppressWarnings(eval(bquote(
        frm(bf(y ~ gp(x, k = .(k))), data = d))))
      bk <- Filter(function(b) identical(b$covstruct, "hsgp"),
                   f$frame$re_blocks)[[1]]
      th <- f$estimates$theta[bk$theta_idx]
      fits[[as.character(k)]] <- list(f = f, th = th)
      dm <- bk$gp_dmax
      L <- 1.25 * 1
      top <- k * pi / (2 * L) / dm
      cat(sprintf(paste0("KS k %3d default start: sigma %.4g sdgp %.4g rho %.4g ",
                         "ll %.3f conv %d msg %s | basis top frequency %.2f\n"),
                  k, sigma(f), exp(th[1]), exp(th[2]) * dm, -f$opt$objective,
                  f$opt$convergence, substr(f$opt$message, 1, 22), top))
    }
    # from k = 60's scales, same data, same k = 30 tape
    st <- fits[["60"]]$th
    for (k in c(15, 30)) {
      f0 <- fits[[as.character(k)]]$f
      bk <- Filter(function(b) identical(b$covstruct, "hsgp"),
                   f0$frame$re_blocks)[[1]]
      tp <- f0$frame$par_template
      tp$theta[bk$theta_idx] <- st
      fn <- frmtmb::build_objective(f0$frame)
      o <- RTMB::MakeADFun(fn, tp, random = "b", map = f0$frame$map, silent = TRUE)
      op <- nlminb(o$par, o$fn, o$gr, control = list(iter.max = 800,
                                                     eval.max = 1200))
      thb <- op$par[names(op$par) == "theta"]
      cat(sprintf(paste0("KS k %3d from k=60's scales: ll %.3f conv %d msg %s, ",
                         "sdgp %.4g rho %.4g | default ll %.3f\n"),
                  k, -op$objective, op$convergence, substr(op$message, 1, 22),
                  exp(thb[1]), exp(thb[2]) * bk$gp_dmax,
                  -f0$opt$objective))
    }
  }
}
