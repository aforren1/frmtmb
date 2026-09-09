# rev-ndt: the eam scale row, timed with the arms INTERLEAVED across
# processes and a control that must report 1.0.
#
# Seed 20260908 for the eam design, seed 99 for the control. One arm per
# process, chosen by REV_ARM:
#
#   new-group      this worktree, ndt_group(s)          the AFTER
#   new-plain      this worktree, no ndt_group()        the BC control
#   old-plain      reflib-r2 (0.6.0), no ndt_group()    the BEFORE
#   new-unbounded  this worktree, max_ndt lifted        Phase 0's arm
#
# The control fit is a gaussian frm() that goes through frmtmb only, and
# frmtmb is the SAME installed package in every arm (this change touches
# frmtmb.eam alone), so its time is the machine and not the model.

arm <- Sys.getenv("REV_ARM", "new-group")
# REV_LIB picks the round: rev-ndt-lib was round one,
# rev-ndt-lib2 the build after the scalar bound went back in the link.
new_lib <- Sys.getenv("REV_LIB",
                      "C:/Users/adf44/source/r/rev-ndt-lib2")
ref_lib <- "C:/Users/adf44/source/r/reflib-r2"
usr_lib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (startsWith(arm, "new")) c(new_lib, ref_lib, usr_lib)
          else c(ref_lib, usr_lib))

suppressMessages({library(frmtmb); library(frmtmb.eam)})

el <- function(expr) {
  gc(FALSE)
  t0 <- Sys.time()
  force(expr)
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
}

# --------------------------------------------------- the control fit
ctrl_time <- el({
  set.seed(99)
  ng <- 30L
  nt <- 400L
  g <- rep(seq_len(ng), each = nt)
  u <- rnorm(ng, 0, 0.7)
  x <- rnorm(ng * nt)
  cd <- data.frame(y = 1 + 0.5 * x + u[g] + rnorm(ng * nt, 0, 1),
                   x = x, g = factor(g))
  cfit <- frm(y ~ x + (1 | g), data = cd, se = TRUE)
})

# ----------------------------------------------------- the eam design
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
                  ndt = tr$ndt * exp(u_nd[s]),
                  bias = 0.5, sv = 0)
d$s <- factor(s)
d$cond <- factor(cond, labels = c("a", "b"))
ndt_true_sub <- tr$ndt * exp(u_nd)

grp <- identical(arm, "new-group")
form <- if (grp) {
  bf(rt | dec(upper) + ndt_group(s) ~ cond + (1 | s), bs ~ 1 + (1 | s),
     ndt ~ 1 + (1 | s), bias = 0.5)
} else {
  bf(rt | dec(upper) ~ cond + (1 | s), bs ~ 1 + (1 | s),
     ndt ~ 1 + (1 | s), bias = 0.5)
}
fam <- if (identical(arm, "new-unbounded")) {
  wiener(max_ndt = 0.45, allow_unreachable = TRUE)
} else wiener()

fit <- NULL
fit_time <- el(fit <- frm(form, family = fam, data = d, se = TRUE))

# ------------------------------------------------------- what it found
dg <- diagnose(fit, quiet = TRUE)
b <- unlist(fixef(fit))
ci <- suppressWarnings(confint(fit))
j <- grep("condb", rownames(ci), fixed = TRUE)[1L]
one <- d[match(levels(d$s), as.character(d$s)), , drop = FALSE]
own <- as.numeric(tapply(d$rt, d$s, min))

nd <- suppressWarnings(
  predict(fit, newdata = d[1L, , drop = FALSE], dpar = "ndt",
          type = "response", re.form = NA, se.fit = TRUE))
frac <- as.numeric(nd$fit[1L])
frac_se <- as.numeric(nd$se.fit[1L])
if (startsWith(arm, "new")) {
  bnd <- frmtmb::single_response(fit)[["family"]][["ndt_bound"]]
  fl <- if (is.null(bnd$floors)) bnd$ub else bnd$floors
  ndt_hat <- frac * mean(fl)
  ndt_se <- frac_se * mean(fl)
  ndt_sub <- suppressWarnings(as.numeric(ndt_time(fit, newdata = one)))
} else {
  ndt_hat <- frac
  ndt_se <- frac_se
  ndt_sub <- suppressWarnings(as.numeric(
    predict(fit, newdata = one, dpar = "ndt", type = "response")))
}
margin_ms <- 1000 * (own - ndt_sub)
vc <- VarCorr(fit)
sds <- vapply(vc, function(m) sqrt(m[1L, 1L]), numeric(1))

out <- c(
  arm = arm,
  ctrl_s = sprintf("%.3f", ctrl_time),
  fit_s = sprintf("%.3f", fit_time),
  logLik = sprintf("%.6f", as.numeric(logLik(fit))),
  conv = as.character(fit$opt$convergence),
  maxgrad = sprintf("%.6g", dg$max_grad),
  pdHess = as.character(isTRUE(dg$pdHess)),
  nbadse = as.character(length(dg$bad_se)),
  npar = as.character(length(fit$opt$par)),
  mu_cond = sprintf("%.6f", unname(b["mu.condb"])),
  mu_cond_lo = sprintf("%.6f", ci[j, 1L]),
  mu_cond_hi = sprintf("%.6f", ci[j, 2L]),
  ndt_pop = sprintf("%.8f", ndt_hat),
  ndt_pop_se = sprintf("%.8f", ndt_se),
  sd_ndt_nat = sprintf("%.8f", sd(ndt_sub)),
  sd_ndt_true = sprintf("%.8f", sd(ndt_true_sub)),
  sd_mu = sprintf("%.6f", sds[1L]),
  sd_bs = sprintf("%.6f", sds[2L]),
  below_own = as.character(sum(margin_ms > 0)),
  n_sub = as.character(length(ndt_sub)),
  min_margin_ms = sprintf("%.4f", min(margin_ms)))

cat("REVROW\t", paste(names(out), out, sep = "=", collapse = "\t"),
    "\n", sep = "")
