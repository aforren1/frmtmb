# Does frmtmb.eam's wiener() fit brms's wiener usage?
#
#   PORT_LIB=<libs> Rscript dev/remeasure-wiener.R
#
# No brms 2.23.0 vignette fits a wiener model: brms_families.Rmd only
# describes it (drift mu, boundary separation bs, non-decision time
# ndt, bias; the response written rt | dec(dec)). So the check is the
# brms spelling of that description, run three ways: brms's stancode()
# (no compile), frm() with frmtmb.eam's wiener(), and an independent
# maximum-likelihood fit through rtdists::ddiffusion() and optim().
# Seed 20261013, n = 400.
lib <- Sys.getenv("PORT_LIB", unset = "")
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})
cat("frmtmb", format(packageVersion("frmtmb")), "frmtmb.eam",
    format(packageVersion("frmtmb.eam")), "from",
    dirname(system.file(package = "frmtmb")), "\n")
set.seed(20261013)
n <- 400
x <- rnorm(n)
mu <- 0.5 + 0.6 * x
d <- do.call(rbind, lapply(seq_len(n), function(i) {
  RWiener::rwiener(1, alpha = 1.4, tau = 0.3, beta = 0.5, delta = mu[i])
}))
dat <- data.frame(rt = d$q, dec = as.integer(d$resp == "upper"), x = x)
f <- bf(rt | dec(dec) ~ x, bs ~ 1, ndt ~ 1, bias ~ 1)

cat("\n## brms 2.23.0 accepts the spelling (stancode, no compile)\n")
sc <- tryCatch(brms::stancode(brms::bf(rt | dec(dec) ~ x, bs ~ 1, ndt ~ 1,
                                       bias ~ 1),
                              family = brms::wiener(), data = dat),
               error = function(e) e)
cat(if (inherits(sc, "error")) paste("ERROR:", conditionMessage(sc)) else
  paste("ok;", sum(grepl("wiener", strsplit(sc, "\n")[[1]])),
        "lines of the Stan code name wiener"), "\n")

cat("\n## frmtmb.eam wiener(), the same formula\n")
fit <- tryCatch(frm(f, family = wiener(), data = dat), error = function(e) e)
if (inherits(fit, "error")) {
  cat("ERROR:", conditionMessage(fit), "\n")
} else {
  print(fixef(fit))
  cat(sprintf("logLik %.6f, convergence %d\n", as.numeric(logLik(fit)),
              fit$opt$convergence))
}

cat("\n## brms's family object handed to frm(), as a brms user would\n")
fb <- tryCatch(frm(f, family = brms::wiener(), data = dat),
               error = function(e) e)
cat(if (inherits(fb, "error")) paste("ERROR:", conditionMessage(fb)) else
  sprintf("ok, logLik %.6f", as.numeric(logLik(fb))), "\n")

cat("\n## independent ML through rtdists::ddiffusion and optim\n")
# RWiener::dwiener() reads delta[1] for every row, so a row-varying
# drift needs rtdists, which agrees with RWiener on scalar arguments.
chk <- c(rtdists::ddiffusion(c(0.6, 0.9), response = c("upper", "lower"),
                             a = 1.4, v = c(0.2, 1.5), t0 = 0.3, z = 0.7),
         mapply(function(q, v, r) RWiener::dwiener(q, alpha = 1.4, tau = 0.3,
                                                   beta = 0.5, delta = v,
                                                   resp = r),
                c(0.6, 0.9), c(0.2, 1.5), c("upper", "lower")))
cat(sprintf("rtdists vs RWiener, two rows: max rel diff %.2g\n",
            max(abs(chk[1:2] / chk[3:4] - 1))))
nll <- function(p) {
  a <- exp(p[3]); tau <- exp(p[4]); b <- plogis(p[5])
  if (tau >= min(dat$rt)) return(1e10)
  -sum(log(rtdists::ddiffusion(dat$rt, response = ifelse(dat$dec == 1,
                                  "upper", "lower"), a = a,
                                v = p[1] + p[2] * dat$x, t0 = tau,
                                z = b * a)))
}
o <- optim(c(0, 0, log(1), log(0.2), 0), nll, method = "BFGS",
           control = list(maxit = 2000, reltol = 1e-12))
o <- optim(o$par, nll, method = "Nelder-Mead",
           control = list(maxit = 5000, reltol = 1e-14))
ref <- c(Intercept = o$par[1], x = o$par[2], bs = exp(o$par[3]),
         ndt = exp(o$par[4]), bias = plogis(o$par[5]))
print(signif(ref, 7))
cat(sprintf("logLik %.6f, convergence %d\n", -o$value, o$convergence))
if (!inherits(fit, "error")) {
  cat(sprintf("logLik difference frmtmb - rtdists ML: %.3g\n",
              as.numeric(logLik(fit)) + o$value))
  fe <- fixef(fit)
  cat(sprintf("drift intercept %.6f vs %.6f, slope %.6f vs %.6f\n",
              fe["Intercept", 1], ref[["Intercept"]], fe["x", 1],
              ref[["x"]]))
}

cat("\n## the extension probe's absent case: an eam export that is not a",
    "family\n")
cat("ddm_simulate exported:", "ddm_simulate" %in%
      getNamespaceExports("frmtmb.eam"), "; accepted as a family:",
    tryCatch(inherits(frmtmb:::as_frmtmb_family(ddm_simulate()),
                      "frmtmb_family"), error = function(e) FALSE), "\n")

cat("\n## the scale of the ndt and bias coefficients against brms's links\n")
# brms's wiener() links: bs log, ndt log, bias logit (brms_families.Rmd
# names the parameters; brms::wiener() carries the links).
bl <- brms::wiener()
cat("brms links: bs", bl$link_bs, " ndt", bl$link_ndt, " bias", bl$link_bias,
    "\n")
if (!inherits(fit, "error")) {
  fe <- fixef(fit)[, 1]
  cat(sprintf("exp(bs_Intercept) %.6f vs reference bs %.6f\n",
              exp(fe[["bs_Intercept"]]), ref[["bs"]]))
  cat(sprintf("exp(ndt_Intercept) %.6f vs reference ndt %.6f\n",
              exp(fe[["ndt_Intercept"]]), ref[["ndt"]]))
  cat(sprintf("plogis(ndt_Intercept) * min(rt) %.6f (min rt %.6f)\n",
              plogis(fe[["ndt_Intercept"]]) * min(dat$rt), min(dat$rt)))
  cat(sprintf("plogis(bias_Intercept) %.6f vs reference bias %.6f\n",
              plogis(fe[["bias_Intercept"]]), ref[["bias"]]))
}
