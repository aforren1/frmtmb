# Lane gpfield: gsigma step 40 (dev/gpstop-ng.R, seed 1): the working
# marginal's other maxima at the fit's optimum, and the objective at the
# better regular maximum's theta with the fit's other values.
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs("lane")
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-rows.R")
ns <- asNamespace("frmtmb")
r <- gf_row("ng:gsigma:step:40:1")
f <- gf_fit(r$fml, data = r$data, family = r$family)$fit
par <- f$opt$par
th <- which(names(par) == "theta")
p <- par; p[th] <- c(log(0.7552), log(0.7996))
cat("objective at the other theta, own coefs", -f$obj$fn(p), "vs", -f$obj$fn(par), "\n")
o <- nlminb(p, f$obj$fn, f$obj$gr)
cat("descent from it", -o$objective, o$convergence, o$par, "\n")
cat("fit pars", par, "\n")
# the working marginal's alternatives
real <- ns$gp_marginal_opt
alts <- NULL
gf_put("gp_marginal_opt", function(s, keep = 3L, distinct = 4L, lrs = 0) {
  m <- .GlobalEnv$real(s, keep, distinct, lrs)
  .GlobalEnv$alts <- m$alts
  m
})
m <- ns$gp_mean_relin(f$obj, f$frame, par)
cat("relin theta", m$theta, "\n")
for (a in alts) {
  p2 <- par; p2[th] <- a
  o2 <- nlminb(p2, f$obj$fn, f$obj$gr)
  cat("alt", a, "objective there", -f$obj$fn(p2), "descent", -o2$objective,
      o2$convergence, "\n")
}
