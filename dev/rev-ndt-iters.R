# rev-ndt: a LOAD-INDEPENDENT count for the speed claim.
#
# The 269 s to 108 s row is two wall clocks on two builds, and a wall
# clock on this machine moves with whatever else is running: over three
# interleaved rounds the control fit, which is the SAME installed
# frmtmb in both arms, ran from 1.44 s to 2.58 s. The optimizer's own
# function and gradient counts do not move with load, and they are the
# actual cause: the global-bound arm does not converge and spends its
# evaluations walking a wall.
#
# REV_ARM = "new-group" | "old-plain". Seed 20260908.

arm <- Sys.getenv("REV_ARM", "new-group")
new_lib <- "C:/Users/adf44/source/r/rev-ndt-lib"
ref_lib <- "C:/Users/adf44/source/r/reflib-r2"
usr_lib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (startsWith(arm, "new")) c(new_lib, ref_lib, usr_lib)
          else c(ref_lib, usr_lib))
suppressMessages({library(frmtmb); library(frmtmb.eam)})

tr <- list(mu0 = 0.4, mu_cond = 0.9, bs = 1.4, ndt = 0.25,
           sd_mu = 0.35, sd_lbs = 0.20, sd_lndt = 0.12)
set.seed(20260908L)
ns <- 30L
nt <- 400L
u_mu <- rnorm(ns, 0, tr$sd_mu)
u_bs <- rnorm(ns, 0, tr$sd_lbs)
u_nd <- rnorm(ns, 0, tr$sd_lndt)
s <- rep(seq_len(ns), each = nt)
cond <- rep(rep(0:1, each = nt / 2L), times = ns)
d <- ddm_simulate(ns * nt,
                  mu = tr$mu0 + tr$mu_cond * cond + u_mu[s],
                  bs = tr$bs * exp(u_bs[s]),
                  ndt = tr$ndt * exp(u_nd[s]), bias = 0.5, sv = 0)
d$s <- factor(s)
d$cond <- factor(cond, labels = c("a", "b"))

form <- if (identical(arm, "new-group")) {
  bf(rt | dec(upper) + ndt_group(s) ~ cond + (1 | s), bs ~ 1 + (1 | s),
     ndt ~ 1 + (1 | s), bias = 0.5)
} else {
  bf(rt | dec(upper) ~ cond + (1 | s), bs ~ 1 + (1 | s),
     ndt ~ 1 + (1 | s), bias = 0.5)
}

t0 <- Sys.time()
fit <- frm(form, family = wiener(), data = d, se = FALSE)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
op <- fit$opt
cat(sprintf(paste0("ITERROW\tarm=%s\tfit_s=%.2f\tconv=%d\t",
                   "logLik=%.6f\titerations=%s\tcounts=%s\n"),
            arm, secs, op$convergence, as.numeric(logLik(fit)),
            paste(op$iterations, collapse = ","),
            paste(names(op$counts), op$counts, sep = "=",
                  collapse = ";")))
