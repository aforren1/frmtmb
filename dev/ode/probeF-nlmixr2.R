# Probe F (opportunistic): cross-check the frmtmb nl-body ODE fit against
# nlmixr2's FOCEi on identical data. The hand-simulated truth is the
# anchor; this only says whether frmtmb lands where the reference PK
# package lands.
source("dev/ode/pk-common.R")
d <- sim_pk()

nd <- rbind(
  data.frame(id = as.integer(unique(d$id)), time = 0, dv = NA,
             amt = PK_TRUTH$dose, evid = 1, cmt = "depot"),
  data.frame(id = as.integer(d$id), time = d$time, dv = d$conc,
             amt = 0, evid = 0, cmt = "cp"))
nd <- nd[order(nd$id, nd$evid == 0, nd$time), ]

library(nlmixr2)
one.cmt <- function() {
  ini({
    tka <- 0
    tke <- log(0.25)
    tv <- log(8)
    eta.ka ~ 0.09
    eta.ke ~ 0.0625
    add.sd <- 0.3
  })
  model({
    ka <- exp(tka + eta.ka)
    ke <- exp(tke + eta.ke)
    v <- exp(tv)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka * depot - ke * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}
tt <- system.time(
  fit <- try(nlmixr2(one.cmt, nd, est = "focei",
                     control = foceiControl(print = 0)), silent = TRUE))
if (inherits(fit, "try-error")) {
  cat("NLMIXR2 FAILED:\n"); cat(as.character(fit), "\n")
} else {
  cat("elapsed:", tt[["elapsed"]], "s\n")
  print(fit)
  cat("\n--- fixed effects ---\n"); print(fixef(fit))
  cat("\n--- omega ---\n"); print(fit$omega)
  cat("\n--- objf ---\n"); print(fit$objDf)
}
