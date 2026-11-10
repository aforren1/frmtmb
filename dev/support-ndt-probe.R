# Lane support, review finding 3: a bound below a shifted support. A
# trunc() lower bound at or below the shift (ndt) is no truncation, so
# the fit must equal the one without it. shifted_lognormal() in the
# reviewer's five forms (seed 402, n 300, y = 2 rexp + 0.05, as
# dev/review-support-trunc.R), then the frmtmb.eam families with an ndt
# that take trunc(): wiener() and rdm() (seed 403), and the slot's taped
# gradient at q below ndt.
#
#   Rscript dev/support-ndt-probe.R   (FRMTMB_PORT_LIB picks the build)
lib <- Sys.getenv("FRMTMB_PORT_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
cmp <- function(label, f0, f1) {
  if (is.character(f1)) {
    cat(sprintf("%-40s ERROR: %s\n", label, substr(f1, 1, 110)))
    return(invisible())
  }
  cat(sprintf("%-40s dlogLik %.3g  max|dbeta| %.3g\n", label,
              as.numeric(logLik(f1)) - as.numeric(logLik(f0)),
              max(abs(f1$opt$par - f0$opt$par))))
}
try_fit <- function(...) {
  tryCatch(suppressWarnings(frm(...)), error = function(e) conditionMessage(e))
}
set.seed(402)
n <- 300
d <- data.frame(y = 2 * stats::rexp(n) + 0.05, lb = -1, lb0 = 0, ub = 1e6)
d$cc <- ifelse(d$y > stats::quantile(d$y, 0.9), 1, 0)
sl <- shifted_lognormal()
f0 <- frm(y ~ 1, family = sl, data = d)
f0c <- frm(y | cens(cc) ~ 1, family = sl, data = d)
cat(sprintf("shifted_lognormal untruncated logLik %.10g\n",
            as.numeric(logLik(f0))))
cmp("shifted_lognormal trunc(lb = -1)",
    f0, try_fit(y | trunc(lb = lb) ~ 1, family = sl, data = d))
cmp("shifted_lognormal trunc(lb = 0)",
    f0, try_fit(y | trunc(lb = lb0) ~ 1, family = sl, data = d))
cmp("shifted_lognormal trunc(lb = -1, ub = 1e6)",
    f0, try_fit(y | trunc(lb = lb, ub = ub) ~ 1, family = sl, data = d))
cmp("shifted_lognormal cens + trunc(lb = -1)",
    f0c, try_fit(y | cens(cc) + trunc(lb = lb) ~ 1, family = sl, data = d))
cmp("shifted_lognormal cens + trunc(lb = 0)",
    f0c, try_fit(y | cens(cc) + trunc(lb = lb0) ~ 1, family = sl, data = d))
# the slot: lS = 0 with a zero gradient below ndt, lF finite
g <- RTMB::MakeTape(function(p) {
  lp <- sl$cdf_log(c(-1, 0), list(mu = p[1], sigma = p[2], ndt = p[3]),
                   list())
  c(lp$lS, lp$lF)
}, c(0, 1, 0.02))
cat("slot at q = -1, 0 (ndt 0.02): value", format(g(c(0, 1, 0.02)),
                                                    digits = 6), "\n")
J <- g$jacobian(c(0, 1, 0.02))
cat("  taped Jacobian rows lS(-1), lS(0), lF(-1), lF(0):",
    format(J, digits = 4), "\n")
cat("  non-finite entries:", sum(!is.finite(J)), "\n")
# frmtmb.eam
set.seed(403)
sw <- ddm_simulate(300, mu = 1, bs = 1.5, ndt = 0.3)
dw <- data.frame(rt = sw$rt, upper = sw$upper,
                 lb = -1, lb0 = 0)
w0 <- frm(bf(rt | dec(upper) ~ 1, bias = 0.5), family = wiener(), data = dw)
cmp("wiener trunc(lb = -1)", w0,
    try_fit(bf(rt | dec(upper) + trunc(lb = lb) ~ 1, bias = 0.5),
            family = wiener(), data = dw))
cmp("wiener trunc(lb = 0)", w0,
    try_fit(bf(rt | dec(upper) + trunc(lb = lb0) ~ 1, bias = 0.5),
            family = wiener(), data = dw))
set.seed(404)
sr <- rdm_simulate(300, v = c(2, 1), ndt = 0.3)
dr <- data.frame(rt = sr$rt, resp = sr$choice, lb = -1, lb0 = 0)
r0 <- try_fit(bf(rt | vint(resp) ~ 1), family = rdm(2), data = dr)
if (is.character(r0)) {
  cat("rdm untruncated ERROR:", substr(r0, 1, 150), "\n")
} else {
  cmp("rdm trunc(lb = -1)", r0,
      try_fit(bf(rt | vint(resp) + trunc(lb = lb) ~ 1), family = rdm(2),
              data = dr))
  cmp("rdm trunc(lb = 0)", r0,
      try_fit(bf(rt | vint(resp) + trunc(lb = lb0) ~ 1), family = rdm(2),
              data = dr))
}
# other bounded supports: royston_parmar() (y > 0) and von_mises()
# (-pi, pi], with a lower bound at or past the edge
suppressMessages(library(frmtmb.survival))
set.seed(405)
ds <- data.frame(y = stats::rweibull(300, 1.5, 2), lb = -1, lb0 = 0)
s0 <- try_fit(y ~ 1, family = royston_parmar(), data = ds)
if (is.character(s0)) {
  cat("royston_parmar untruncated ERROR:", substr(s0, 1, 150), "\n")
} else {
  cmp("royston_parmar trunc(lb = -1)", s0,
      try_fit(y | trunc(lb = lb) ~ 1, family = royston_parmar(), data = ds))
  cmp("royston_parmar trunc(lb = 0)", s0,
      try_fit(y | trunc(lb = lb0) ~ 1, family = royston_parmar(), data = ds))
}
set.seed(406)
dv <- data.frame(y = atan2(sin(stats::rnorm(300, 0.5, 0.8)),
                           cos(stats::rnorm(300, 0.5, 0.8))), lb = -4,
                 lbp = -pi)
v0 <- frm(y ~ 1, family = von_mises(), data = dv)
cmp("von_mises trunc(lb = -4)", v0,
    try_fit(y | trunc(lb = lb) ~ 1, family = von_mises(), data = dv))
cmp("von_mises trunc(lb = -pi)", v0,
    try_fit(y | trunc(lb = lbp) ~ 1, family = von_mises(), data = dv))
